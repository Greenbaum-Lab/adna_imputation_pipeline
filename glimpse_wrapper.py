# IMPORTS
import os
import re
import subprocess
from concordance_plot import plot
from utils import (validate_chromosome_notation_in_bam, CHR_NOTATION_REGEX,
                   create_bcf_from_bam)

# CONSTS
CHR_NOTATION_REGEX = '^chr([1-9]|1[0-9]|2[012]|X|Y|M|MT)$'
MAP = 'map/genetic_maps.b37/{chromosome}.b37.gmap.gz'
WD_REGEX = '(\w+_subsample_\d+)'
CHUNK_OUTPUT = 'chunks.txt'
PHASE_OUTPUT_DIR = 'GLIMPSE_impute'
LST_FILE_NAME = 'list.txt'
LIGATE_OUTPUT = 'ligated.bcf'


class GlimpseException(Exception):
    ERROR_MESSAGE = ("GLIMPSE {} DID NOT RUN SUCCESSFULLY, "
                     "NO OUTPUT FILE WAS CREATED.")

    def __init__(self, stage_name):
        message = GlimpseException.ERROR_MESSAGE.format(stage_name)
        super().__init__(message)


class GlimpseWrapper(object):
    def __init__(self, glimpse_dir, wd, input_bam, chromosome, reference,
                 truth_fasta, reference_sites=None, truth_bam=None, truth_bcf=None,
                 allele_freq_file=None):
        self.glimpse_dir = glimpse_dir
        self.wd = wd
        self.input_bam = input_bam
        # TODO maybe fix notation instead of throwing an error?
        if not self.validate_chromosome_notation(chromosome):
            raise ValueError('Chromosome must be of notation chr<1-23/X/Y/MT/M>')
        self.chromosome = chromosome
        self.reference = reference
        if not os.path.isfile('.'.join([reference, 'csi'])):
            subprocess.run(['bcftools', 'index', '-f', self.reference])
        self.map_file = os.path.join(wd, MAP.format(chromosome=chromosome))
        self.truth_fasta = truth_fasta
        self.reference_sites = reference_sites
        self.truth_bam = truth_bam
        self.truth_bcf = truth_bcf
        self.allele_freq_file = allele_freq_file

        if not os.path.isdir(os.path.join(wd, 'output')):
            subprocess.run(['mkdir', os.path.join(wd, 'output')])
        self.output_dir = os.path.join(wd, 'output',
                                       os.path.splitext(os.path.basename(input_bam))[0])
        subprocess.run(['mkdir', self.output_dir])
        with open(f'{self.output_dir}/params.txt', 'w') as f:
            f.write(f'input_bam: {self.input_bam}\n'
                    f'chromosome: {self.chromosome}\n'
                    f'reference: {self.reference}\n'
                    f'map_file: {self.map_file}\n'
                    f'truth_fasta: {self.truth_fasta}\n'
                    f'reference_sites: {self.reference_sites}\n'
                    f'truth_bam: {self.truth_bam}\n'
                    f'truth_bcf: {self.truth_bcf}\n'
                    f'allele_freq_file: {self.allele_freq_file}\n')

    @staticmethod
    def validate_chromosome_notation(chromosome):
        """
        This method makes sure the given chromosome matches hg38 notation - chrXX.

        :param chromosome: Given chromosome
        :return: True - if the chromosome matches the notation
                 False - otherwise
        """
        return re.compile(CHR_NOTATION_REGEX).match(chromosome) is not None

    def extract_sites_from_reference(self):
        """
        This method extracts sites from the reference file for future steps.
        """
        output_path = (self.reference[:self.reference.rfind('.')] +
                       '.sites.vcf.gz')
        command = subprocess.run(['bcftools', 'view', '-G', '-Oz',  '-o',
                                  output_path, self.reference],
                                 capture_output=True, text=True)
        command = subprocess.run(['bcftools', 'index', '-f', output_path],
                                 capture_output=True, text=True)
        self.reference_sites = output_path

    def run_chunk(self):
        """
        This method runs chunk step of glimpse, i.e. splitting the extracted
        sites from the reference into chunks.
        """
        if (not self.reference_sites and
                not os.path.isfile(self.reference[:self.reference.rfind('.')] +
                                   '.sites.vcf.gz')):
            self.extract_sites_from_reference()
        elif not self.reference_sites:
            self.reference_sites = (self.reference[:self.reference.rfind('.')] +
                                    '.sites.vcf.gz')

        chunk_stage_path = os.path.join(self.glimpse_dir,
                                        'chunk/bin/GLIMPSE2_chunk')
        output_path = os.path.join(os.path.dirname(self.reference), CHUNK_OUTPUT)
        if os.path.isfile(output_path):
            print('-------------------CHUNK------------------\n' +
                  'chunk output file already exists')
            return

        command = subprocess.Popen([chunk_stage_path, '--input',
                                    self.reference_sites, '--region',
                                    self.chromosome, '--output', output_path,
                                    '--map', self.map_file, '--sequential'],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return_code = command.wait()
        stdout, stderr = command.communicate()
        if return_code:
            raise Exception(stderr)

        print('-------------------CHUNK------------------\n', stdout)
        if not os.path.isfile(output_path):
            raise GlimpseException("CHUNK")

    def run_split_reference(self):
        """
        This method runs split_reference step of glimpse, i.e. converting the
        reference panel into GLIMPSE2’s binary file format.
        """
        split_reference_stage_path = os.path.join(self.glimpse_dir,
                                                  'split_reference/bin/'
                                                  'GLIMPSE2_split_reference')

        reference_basename = os.path.basename(self.reference)[:os.path.basename(self.reference).rfind('.')]
        split_reference_output_path = os.path.join(os.path.dirname(self.reference),
                                                   f'{reference_basename}_split')
        subprocess.run(['mkdir', '-p', split_reference_output_path])

        with open(os.path.join(os.path.dirname(self.reference), CHUNK_OUTPUT), "r") as file:
            for line in file:
                line = line.strip()
                if line:
                    parts = line.split()
                    irg_value = parts[2]
                    org_value = parts[3]

                    output = os.path.join(split_reference_output_path, f'{reference_basename}')

                    split_range = irg_value.split(':')[1].split('-')
                    if os.path.isfile(f'{output}_{self.chromosome}_'
                                      f'{split_range[0]}_{split_range[1]}.bin'):
                        print('-------------------SPLIT_REF------------------\n'
                              f'bin file for range {split_range[0]}-{split_range[1]}'
                              f' already exists')
                        continue
                    # Construct the command to run
                    command = [
                        split_reference_stage_path,
                        "--reference", self.reference,
                        "--map", self.map_file,
                        "--input-region", irg_value,
                        "--output-region", org_value,
                        "--output", output
                    ]

                    processed = subprocess.run(command, capture_output=True, text=True)
                    # TODO check if the command succeeded
                    print("-------------------SPLIT_REF------------------\n" + processed.stdout)

    def filter_only_relevant_chromosome_bam(self, bam_file):
        """
        This method filters only the wanted chromosome from the BAM file,
        in order to make the next steps run faster.
        """
        # TODO perhaps this check lies because the header stays the same after filtering
        print(f"Running command: " + "samtools view -H " + bam_file +
                                     " | grep '^@SQ' | awk '{print $2}' | cut -f 2 -d ':'")
        chromosomes = subprocess.run("samtools view -H " + bam_file +
                                     " | grep '^@SQ' | awk '{print $2}' | cut -f 2 -d ':'",
                                     shell=True, capture_output=True, text=True).stdout.splitlines()
        # This check assumes the chromosomes notation matches self.chromosome
        if len(chromosomes) == 1 and chromosomes[0] == self.chromosome:
            return bam_file
        original_bam_name = os.path.basename(bam_file)
        new_bam_name = (original_bam_name[:original_bam_name.rfind('.')] +
                        f'_{self.chromosome}.bam')
        new_bam_path = os.path.join(self.output_dir, new_bam_name)
        if os.path.isfile(new_bam_path):
            print('Filtered BAM file already exists')
        else:
            print(f'Running command: samtools view -bh {bam_file} {self.chromosome} > {new_bam_path}')
            command = subprocess.run(f'samtools view -bh {bam_file} '
                                     f'{self.chromosome} > {new_bam_path}',
                                     capture_output=True, text=True, shell=True)
            print(f'Running command: samtools index {new_bam_path}')
            command = subprocess.run(['samtools', 'index', new_bam_path],
                                     capture_output=True, text=True)
        return new_bam_path

    def run_phase(self):
        """
        This method runs the phase step of glimpse.
        """
        # Making sure chromosome notation matches between reference and input
        # BAM file
        # TODO change order of those functions (requires changes in filter function)
        self.input_bam = validate_chromosome_notation_in_bam(self.input_bam,
                                                             self.output_dir,
                                                             self.chromosome)
        # TODO I need to rethink how I check if this function is needed for a given file
        # Commented out this part since I assume the BAM files will be filtered
        # before running the pipeline.
        # self.input_bam = self.filter_only_relevant_chromosome_bam(self.input_bam)

        subprocess.run(['mkdir', os.path.join(self.output_dir, PHASE_OUTPUT_DIR)])
        reference_basename = os.path.basename(self.reference)[:os.path.basename(self.reference).rfind('.')]
        reference_files = os.listdir(os.path.join(os.path.dirname(self.reference),
                                                  f'{reference_basename}_split'))

        phase_stage_path = os.path.join(self.glimpse_dir,
                                        'phase/bin/GLIMPSE2_phase')

        for reference_file in reference_files:
            output_file = 'imputed_' + reference_file[:reference_file.rfind('.')] + '.bcf'
            output_path = os.path.join(self.output_dir, PHASE_OUTPUT_DIR,
                                       output_file)
            if os.path.isfile(output_path):
                print('-------------------PHASE------------------\n' +
                      f'phased file {output_file} already exists')
                continue

            command = [
                phase_stage_path,
                '--bam-file', self.input_bam,
                '--reference',
                os.path.join(os.path.dirname(self.reference),
                             f'{reference_basename}_split',
                             reference_file),
                '--output', output_path
            ]

            processed = subprocess.run(command, capture_output=True, text=True)
            # TODO make sure the stage succeeded
            print("-------------------PHASE------------------\n" + processed.stdout)

    def run_ligate(self):
        """
        This method runs the ligate step of glimpse.
        """
        command = subprocess.run(f'ls -1v {os.path.join(self.output_dir, PHASE_OUTPUT_DIR)}'
                                 f'/imputed_*.bcf > {os.path.join(self.output_dir, LST_FILE_NAME)}',
                                 shell=True, capture_output=True, text=True)
        ligate_stage_path = os.path.join(self.glimpse_dir,
                                         'ligate/bin/GLIMPSE2_ligate')
        ligated_file_name = (f'{os.path.basename(self.input_bam)[:os.path.basename(self.input_bam).rfind(".")]}'
                             f'_{LIGATE_OUTPUT}')
        if os.path.isfile(os.path.join(self.output_dir, ligated_file_name)):
            print('-------------------LIGATE------------------\n' +
                  'ligated file already exists')
            return
        command = subprocess.run([ligate_stage_path, '--input',
                                  os.path.join(self.output_dir, LST_FILE_NAME),
                                  '--output', os.path.join(self.output_dir, LIGATE_OUTPUT)],
                                 capture_output=True, text=True)
        print("-------------------LIGATE------------------\n" + command.stdout)

    def create_bcf_from_bam(self):
        """
        This method creates a BCF file deriving from the BAM file in order to
        run concordance step.
        """
        bcf_name = self.truth_bam[:self.truth_bam.rfind('.')].replace('bam', 'bcf') + '.bcf'
        create_bcf_from_bam(input_bam=self.truth_bam,
                            chromosome=self.chromosome,
                            truth_fasta=self.truth_fasta,
                            reference_sites=self.reference_sites,
                            output_path=bcf_name)
        self.truth_bcf = bcf_name

    def run_concordance(self):
        """
        This method runs the concordance step of glimpse, i.e. extracting
        information about the imputation's accuracy. To be run only if a
        truth_bam file was provided.
        """
        if not self.allele_freq_file:
            self.allele_freq_file = self.reference_sites

        if not self.truth_bcf:
            print(f'Current truth_bam: {self.truth_bam}')
            self.truth_bam = validate_chromosome_notation_in_bam(self.truth_bam,
                                                                 self.output_dir,
                                                                 self.chromosome)
            print(f'Current truth_bam: {self.truth_bam}')
            # TODO I need to rethink how I check if this function is needed for a given file
            self.truth_bam = self.filter_only_relevant_chromosome_bam(self.truth_bam)
            print(f'Current truth_bam: {self.truth_bam}')
            self.create_bcf_from_bam()

        # Validate sample names match between BCF files
        header = subprocess.run(['bcftools', 'view', '-h',
                                 os.path.join(self.output_dir, LIGATE_OUTPUT)],
                                capture_output=True, text=True).stdout
        truth_sample_name = subprocess.run(['bcftools', 'query', '-l',
                                           self.truth_bcf], capture_output=True,
                                           text=True).stdout
        header = header[:header.rfind('\t')] + '\t' + truth_sample_name
        with open(os.path.join(self.output_dir, 'fixed_header.txt'), 'w') as f:
            f.write(header)
        renamed_ligated_output = os.path.join(self.output_dir,
                                              'renamed_' + LIGATE_OUTPUT)
        command = subprocess.run(['bcftools', 'reheader', '-h',
                                  os.path.join(self.output_dir, 'fixed_header.txt'),
                                  os.path.join(self.output_dir, LIGATE_OUTPUT),
                                  '-o', renamed_ligated_output],
                                 capture_output=True, text=True)
        command = subprocess.run(['bcftools', 'index', '-f',
                                  renamed_ligated_output],
                                 capture_output=True, text=True)

        # Creating concordance.lst file
        subprocess.run(f'echo {self.chromosome} {self.allele_freq_file} '
                       f'{self.truth_bcf} {renamed_ligated_output} > '
                       f'{self.output_dir}/concordance.lst', shell=True)

        concordance_stage_path = os.path.join(self.glimpse_dir,
                                              'concordance/bin/GLIMPSE2_concordance')

        subprocess.run(['mkdir', os.path.join(self.output_dir, 'GLIMPSE_concordance')])
        if os.listdir(os.path.join(self.output_dir, 'GLIMPSE_concordance')):
            print("-------------------CONCORDANCE------------------\n" +
                  'concordance directory exists and is not empty')
            return

        command = subprocess.run([concordance_stage_path, '--input',
                                  os.path.join(self.output_dir, 'concordance.lst'),
                                  '--min-val-dp', '8', '--output',
                                  os.path.join(self.output_dir, 'GLIMPSE_concordance/output'),
                                  '--min-val-gl', '0.9999', '--bins', '0.00000',
                                  '0.00100', '0.00200', '0.00500', '0.01000',
                                  '0.05000', '0.10000', '0.20000', '0.50000',
                                  '--af-tag', 'AF_nfe', '--thread', '4'],
                                 capture_output=True, text=True)
        # TODO make sure the step succeeded
        print("-------------------CONCORDANCE------------------\n" + command.stdout)

        plot(data=os.path.join(self.output_dir, 'GLIMPSE_concordance/output.rsquare.grp.txt.gz'),
             output=os.path.join(self.output_dir, 'accplot.png'))

    def run(self):
        self.run_chunk()
        self.run_split_reference()
        self.run_phase()
        self.run_ligate()
        if self.truth_bam or self.truth_bcf:
            self.run_concordance()

    def get_output_path(self):
        return os.path.join(self.output_dir, LIGATE_OUTPUT)

    def delete_temp_files(self):
        """
        This function deletes all temp files created by the GLIMPSE pipeline.
        That is, only leaving the output of ligate stage and changing it's
        name from renamed.bcf to the sample's name.
        """
        # Remove phase temp files
        subprocess.run(f'rm -r {os.path.join(self.output_dir, PHASE_OUTPUT_DIR)}',
                       shell=True)
        # Remove BAM with changed chromosome notation
        # KEITH: commented out the rm command for the moment to use the chr notation files afterwards for calling
        new_bam_name = (os.path.basename(self.input_bam)[:os.path.basename(self.input_bam).rfind('.')]
                        + '_chr_notation.bam')
        # subprocess.run(['rm', os.path.join(self.output_dir, new_bam_name)])
        subprocess.run(['rm', os.path.join(self.output_dir, 'header_chr_notation.txt')])
        # Remove ligate stage temp file
        subprocess.run(['rm', os.path.join(self.output_dir, LST_FILE_NAME)])


if __name__ == '__main__':
    samples = {'VLASA7': '/home/lab-heavy2/adna/bam/shotgun/VLASA7/VLASA7.bam',
               'Klein7': '/home/lab-heavy2/adna/bam/shotgun/Klein7/Klein7.bam',
               'I0676': '/home/lab-heavy2/adna/bam/1240K/I0676/I0676.bam'}
    for sample_name, sample_path in samples.items():
        print(f'WORKING ON SAMPLE {sample_name}')
        GlimpseWrapper(glimpse_dir='/home/lab-heavy2/Documents/yuvalt/glimpse',
                       wd='/home/lab-heavy2/adna',
                       input_bam=sample_path,
                       chromosome='chr6',
                       reference='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.bcf',
                       reference_sites='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.gz',
                       truth_fasta='/home/lab-heavy2/adna/reference/hg19.fa').run()
    # for level in ['0001', '00025', '0005', '00075', '001', '002']:
    # # for level in ['1']:
    #     print(f'WORKING ON SUBSAMPLE LEVEL OF {level}')
    #     GlimpseWrapper(glimpse_dir='/home/lab-heavy2/Documents/yuvalt/glimpse',
    #                    wd='/home/lab-heavy2/adna',
    #                    # input_bam=f'/home/lab-heavy2/adna/bam/shotgun/Yana_old/Yana_old_chr_notation_chr22.bam',
    #                    input_bam=f'/home/lab-heavy2/adna/bam/1240K/I5241/I5241_chr_notation_chr6_subsampled_{level}x.bam',
    #                    # input_bam=f'/home/lab-heavy2/adna/bam/shotgun/Yana_old/Yana_old_chr_notation_chr22_subsampled_{level}_percent.bam',
    #                    chromosome='chr6',
    #                    reference='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.bcf',
    #                    reference_sites='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.gz',
    #                    # map_file='/home/lab-heavy2/Documents/yuvalt/glimpse/maps/genetic_maps.b37/chr22.b37.gmap.gz',
    #                    truth_fasta='/home/lab-heavy2/adna/reference/hg19.fa',
    #                    # allele_freq_file='/home/lab-heavy2/Documents/yuvalt/glimpse/'
    #                    #                      'tutorial/GLIMPSE_validation/gnomad.genomes'
    #                    #                      '.r3.0.sites.chr22.isec.bcf',
    #                    # # reference_sites='/home/lab-heavy2/adna/reference/1K/1K.chr22/1K.chr22.sites.vcf.gz',
    #                    truth_bam='/home/lab-heavy2/adna/bam/shotgun/VLASA32/VLASA32_chr_notation_chr6.bam',
    #                    truth_bcf='/home/lab-heavy2/adna/bcf/1240K/I5241/I5241_chr_notation_chr6.bcf').run()


