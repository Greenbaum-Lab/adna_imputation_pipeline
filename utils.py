import re
from pathlib import Path
import subprocess
import random
import pysam
import json
import glob
import os
from concordance_plot import plot

COMPRESSED_VCF_EXT = '.vcf.gz'
_1240K_SITES = '/home/lab-heavy2/adna/reference/1240K_positions.lst'
INC_1240K_NAME_ADDITION = 'only_1240K_sites'
EXC_1240K_NAME_ADDITION = 'exclude_1240K_sites'
AREA_ALL = 'entire_chr6'
AREA_MHC = 'MHC_region'
AREA_EXC_MHC = 'chr6_excluding_MHC_region'
CHR_NOTATION_REGEX = '^chr([1-9]|1[0-9]|2[012]|X|Y|M|MT)$'
CHR_NOTATION_REPLACEMENTS = {
                'SN:1': 'SN:chr1',
                'SN:2': 'SN:chr2',
                'SN:3': 'SN:chr3',
                'SN:4': 'SN:chr4',
                'SN:5': 'SN:chr5',
                'SN:6': 'SN:chr6',
                'SN:7': 'SN:chr7',
                'SN:8': 'SN:chr8',
                'SN:9': 'SN:chr9',
                'SN:10': 'SN:chr10',
                'SN:11': 'SN:chr11',
                'SN:12': 'SN:chr12',
                'SN:13': 'SN:chr13',
                'SN:14': 'SN:chr14',
                'SN:15': 'SN:chr15',
                'SN:16': 'SN:chr16',
                'SN:17': 'SN:chr17',
                'SN:18': 'SN:chr18',
                'SN:19': 'SN:chr19',
                'SN:20': 'SN:chr20',
                'SN:21': 'SN:chr21',
                'SN:22': 'SN:chr22',
                'SN:X': 'SN:chrX',
                'SN:Y': 'SN:chrY',
                'SN:M': 'SN:chrM'
            }


def split_chr6_vcf(original_vcf_path, mhc_region_file_path, chr6_exclude_mhc_file_path):
    """
    This function splits a given chr6 compressed VCF file to 2 files - one containing
    only the MHC region, and one containing anything but the MHC region.
    The files will be created in the original file's directory.

    :param original_vcf_path: file to split (compressed VCF)
    :param mhc_region_file_path: path of MHC region compressed VCF file
    :param chr6_exclude_mhc_file_path: path of chr6 excluding MHC region compressed VCF name
    """
    # Input validation
    if (not original_vcf_path.endswith(COMPRESSED_VCF_EXT) or
            not mhc_region_file_path.endswith(COMPRESSED_VCF_EXT) or
            not chr6_exclude_mhc_file_path.endswith(COMPRESSED_VCF_EXT)):
        raise ValueError("Input and output files must have .vcf.gz extensions.")

    # Creating MHC region vcf
    if os.path.isfile(mhc_region_file_path):
        print(f"File {mhc_region_file_path} already exists")
    else:
        subprocess.run(['bcftools', 'view', '-O', 'z', '-r', 'chr6:28477797-33448354',
                        original_vcf_path, '-o', mhc_region_file_path],
                       capture_output=True, text=True)
        subprocess.run(['bcftools', 'index', '-f', '-t', mhc_region_file_path],
                       capture_output=True, text=True)

    # Creating chr6 excluding MHC region vcf
    if os.path.isfile(chr6_exclude_mhc_file_path):
        print(f"File {chr6_exclude_mhc_file_path} already exists")
    else:
        chr6_excluding_mhc_part1_path = chr6_exclude_mhc_file_path.replace(COMPRESSED_VCF_EXT,
                                                                           '_part1' + COMPRESSED_VCF_EXT)
        subprocess.run(['bcftools', 'view', '-O', 'z', '-r', 'chr6:1-28477796',
                        original_vcf_path, '-o', chr6_excluding_mhc_part1_path],
                       capture_output=True, text=True)
        subprocess.run(['bcftools', 'index', '-f', '-t', chr6_excluding_mhc_part1_path],
                       capture_output=True, text=True)

        chr6_excluding_mhc_part2_path = chr6_exclude_mhc_file_path.replace(COMPRESSED_VCF_EXT,
                                                                           '_part2' + COMPRESSED_VCF_EXT)
        subprocess.run(['bcftools', 'view', '-O', 'z', '-r', 'chr6:33448355-',
                        original_vcf_path, '-o', chr6_excluding_mhc_part2_path],
                       capture_output=True, text=True)
        subprocess.run(['bcftools', 'index', '-f', '-t', chr6_excluding_mhc_part2_path],
                       capture_output=True, text=True)

        subprocess.run(['bcftools', 'concat', chr6_excluding_mhc_part1_path,
                        chr6_excluding_mhc_part2_path, '-O', 'z', '-o',
                        chr6_exclude_mhc_file_path],
                       capture_output=True, text=True)
        subprocess.run(['bcftools', 'index', '-f', '-t', chr6_exclude_mhc_file_path],
                       capture_output=True, text=True)

        # Cleanup
        subprocess.run([f'rm {chr6_excluding_mhc_part1_path}*'],
                       capture_output=True, text=True, shell=True)
        subprocess.run([f'rm {chr6_excluding_mhc_part2_path}*'],
                       capture_output=True, text=True, shell=True)


def run_only_concordance(input_bcf, chromosome, allele_freq_file, truth_bcf):
    """
    This function runs only the concordance step (including the creation of temp
    files), without the need to run the entire GLIMPSE pipeline.

    :param input_bcf: to run concordance on
    :param chromosome:
    :param allele_freq_file:
    :param truth_bcf:
    """
    output_dir = os.path.dirname(input_bcf)
    file_name = os.path.basename(input_bcf)[:-4]
    # Creating concordance.lst file
    subprocess.run(f'echo {chromosome} {allele_freq_file} '
                   f'{truth_bcf} {input_bcf} > '
                   f'{output_dir}/concordance_{file_name}.lst', shell=True)

    subprocess.run(['mkdir', os.path.join(output_dir, f'GLIMPSE_concordance_{file_name}')])
    if os.listdir(os.path.join(output_dir, f'GLIMPSE_concordance_{file_name}')):
        print("-------------------CONCORDANCE------------------\n" +
              'concordance directory exists and is not empty')
        return

    command = subprocess.run(['/home/lab-heavy2/Documents/yuvalt/glimpse/concordance/bin/GLIMPSE2_concordance',
                              '--input', os.path.join(output_dir, f'concordance_{file_name}.lst'),
                              '--min-val-dp', '8', '--output',
                              os.path.join(output_dir, f'GLIMPSE_concordance_{file_name}/output'),
                              '--min-val-gl', '0.9999', '--bins', '0.00000',
                              '0.00100', '0.00200', '0.00500', '0.01000',
                              '0.05000', '0.10000', '0.20000', '0.50000',
                              '--af-tag', 'AF_nfe', '--thread', '4'],
                             capture_output=True, text=True)
    # TODO make sure the step succeeded
    print("-------------------CONCORDANCE------------------\n" + command.stdout)


def downsample(input_bam, depths, original_average_depth):
    """
    This function creates downsampled files for a given input BAM file,
    in the BAM file's directory.

    :param input_bam: sample to downsample
    :param depths: wanted average depths of downsampled BAM files
    :param original_average_depth: of input_bam
    """
    for depth in depths:
        seed = random.randint(0, 2 ** 31 - 1)
        cmd = (f'samtools view -b -s {seed + round(depth / original_average_depth, 4)}'
               f' {input_bam} > {input_bam[:-4]}_subsampled_{str(depth).replace(".", "")}x.bam')
        subprocess.run(cmd, shell=True)
        cmd = f'samtools index {input_bam[:-4]}_subsampled_{str(depth).replace(".", "")}x.bam'
        subprocess.run(cmd, shell=True)


def create_compressed_vcf_from_bam(input_bam, chromosome, truth_fasta,
                                   reference_sites, output_path, tsv):
    """
    This function creates a compressed VCF file out of a given BAM file.
    """
    if not output_path.endswith('.vcf.gz'):
        raise ValueError("Output path must end with '.vcf.gz' in order to match "
                         "the output's format.")
    print(f"Running command: bcftools mpileup -f {truth_fasta} -I -E -a "
          f"'FORMAT/DP' -T {reference_sites} -r {chromosome} {input_bam} -Ou | "
          f"bcftools call -Aim -C alleles -T {tsv} -Oz -o {output_path}")
    command = subprocess.run(f"bcftools mpileup -f {truth_fasta} -I -E"
                             f" -a 'FORMAT/DP' -T {reference_sites} -r {chromosome}"
                             f" {input_bam} -Ou | bcftools call -Aim "
                             f"-C alleles -T {tsv} -Oz -o {output_path}",
                             shell=True, capture_output=True, text=True)
    print(f"Running command: bcftools index {output_path}")
    command = subprocess.run(['bcftools', 'index', output_path],
                             capture_output=True, text=True)


def create_bcf_from_bam(input_bam, chromosome, truth_fasta, reference_sites,
                        output_path, tsv=None):
    """
    This function creates a compressed VCF file out of a given BAM file.
    """
    if not output_path.endswith('.bcf'):
        raise ValueError("Output path must end with '.bcf' in order to match "
                         "the output's format.")

    if os.path.isfile(output_path):
        print('BCF file already exists')

    else:
        # Create TSV from VCF
        if not tsv:
            tsv = reference_sites[:reference_sites.rfind('.')] + '.tsv.gz'
            if not os.path.isfile(tsv):
                print(f"Running command: bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n'"
                      f" {reference_sites} | bgzip -c > {tsv} && tabix -s1 -b2 -e2 {tsv}")
            command = subprocess.run(f"bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' "
                                     f"{reference_sites} | bgzip -c > {tsv}"
                                     f" && tabix -s1 -b2 -e2 {tsv}", shell=True,
                                     capture_output=True, text=True)
        print(f"Running command: bcftools mpileup -f {truth_fasta} -I -E -a "
              f"'FORMAT/DP' -T {reference_sites} -r {chromosome} {input_bam}"
              f" -Ou | bcftools call -Aim -C alleles -T {tsv} -Ob -o {output_path}")
        command = subprocess.run(f"bcftools mpileup -f {truth_fasta} -I -E"
                                 f" -a 'FORMAT/DP' -T {reference_sites} -r {chromosome}"
                                 f" {input_bam} -Ou | bcftools call -Aim "
                                 f"-C alleles -T {tsv} -Ob -o {output_path}",
                                 shell=True, capture_output=True, text=True)
        print(f"Running command: bcftools index {output_path}")
        command = subprocess.run(['bcftools', 'index', output_path],
                                 capture_output=True, text=True)


def create_compressed_vcf_from_bcf(input_bcf, output_path):
    """
    This function creates a compressed VCF file out of a given BCF file.

    :param input_bcf: File to convert.
    :param output_path: Path to save converted file in.
    """
    if not output_path.endswith('.vcf.gz'):
        raise ValueError("Output path must end with '.vcf.gz' in order to match "
                         "the output's format.")
    print(f"Running command: bcftools view {input_bcf} -Oz -o {output_path}")
    command = subprocess.run(['bcftools', 'view', input_bcf, '-Oz', '-o', output_path],
                             capture_output=True, text=True)
    print(f"Running command: bcftools index {output_path}")
    command = subprocess.run(['bcftools', 'index', output_path],
                             capture_output=True, text=True)


def filter_1240k_sites_from_vcf_files_in_dir(root_dir, files_name, include):
    """
    This function fetches all files under the given dir recursively and filters
    1240K sites from each of them. Filtered output file will be saved in the
    same location of the input file.

    :param root_dir: Directory to fetch files from.
    :param files_name: Files to fetch.
    :param include: whether to include or exclude 1240k sites
    """
    pattern = f'{root_dir}/**/{files_name}'
    matching_files = glob.glob(pattern, recursive=True)

    for file_path in matching_files:
        name_addition = INC_1240K_NAME_ADDITION if include \
            else EXC_1240K_NAME_ADDITION
        output_path = (file_path[:file_path.find(COMPRESSED_VCF_EXT)] + '_' +
                       name_addition + COMPRESSED_VCF_EXT)
        if os.path.isfile(output_path):
            print(f"File {output_path} already exists")
        else:
            command = subprocess.run(['bcftools', 'view', '-T',
                                      f'{"^" if not include else ""}{_1240K_SITES}',
                                      '-Oz', '-o', output_path, file_path],
                                     capture_output=True, text=True)
            stdout = command.stdout
            command = subprocess.run(['tabix', '-p', 'vcf', output_path],
                                     capture_output=True, text=True)


def analyze_snp_depths_in_vcf(vcf_file, min_dp):
    """
    This function counts the number of sites with distinguish genotype,
    the number of sites without DP information,and the number of sites
    with DP value greater or equal to the given dp, in the given VCF file.

    :param vcf_path: File to count sites of
    :param min_dp: minimum depth to include
    :return: the counts in the order of the function's description.
    """
    vcf = pysam.VariantFile(vcf_file)
    total_count = 0
    without_dp_count = 0
    dp_greater_than_min_dp = 0

    for record in vcf:
        sample = record.samples[0]  # Access the first sample
        if sample.alleles:  # Check if alleles information is available
            gt = sample['GT']
            if gt and None not in gt:  # Check for valid genotype data
                total_count += 1
                dp = sample.get('DP')
                if dp and dp >= min_dp:
                    dp_greater_than_min_dp += 1
                elif not dp:
                    without_dp_count += 1

    vcf.close()
    return total_count, without_dp_count, dp_greater_than_min_dp


def validate_chromosome_notation_in_bam(bam_file, output_dir, chromosome):
    """
    This method validates that the chromosome's notation in the BAM file
    matches the reference.

    :return: new bam name
    """
    # Check if there's need to reheader
    print(f'Running command: samtools view -H {bam_file} | grep "@SQ"')
    output = subprocess.run(f'samtools view -H {bam_file} | grep "@SQ"',
                            shell=True, capture_output=True, text=True).stdout

    # TODO right now I reheader if there's any chromosome with the wrong notation, not the relevant one.
    if chromosome not in output:
        # Convert chromosomes names
        new_bam_name = (os.path.basename(bam_file)[:os.path.basename(bam_file).rfind('.')]
                        + '_chr_notation.bam')
        new_bam_path = os.path.join(output_dir, new_bam_name)
        if os.path.isfile(new_bam_path):
            print('BAM file with the correct chromosome notation already exists')
        else:
            print(f'Running command: samtools view -H {bam_file}')
            bam_header = subprocess.run(['samtools', 'view', '-H', bam_file],
                                        capture_output=True, text=True).stdout
            # Apply the replacements
            for old_sn, new_sn in CHR_NOTATION_REPLACEMENTS.items():
                bam_header = bam_header.replace(old_sn, new_sn)

            header_file = os.path.join(output_dir, 'header_chr_notation.txt')
            with open(header_file, 'w') as f:
                f.write(bam_header)
            print(f'Running command: samtools reheader {header_file} {bam_file} > {new_bam_path}')
            command = subprocess.run(f'samtools reheader {header_file} {bam_file}'
                                     f' > {new_bam_path}',
                                     shell=True, capture_output=True, text=True)
            print(f'Running command: samtools index {new_bam_path}')
            command = subprocess.run(['samtools', 'index', new_bam_path],
                                     capture_output=True, text=True)
        return new_bam_path
    return bam_file


if __name__ == '__main__':
    # run_only_concordance('/home/lab-heavy2/adna/output/simulated_SC_sorted/simulated_SC_MHC_region.bcf',
    #                      'chr6',
    #                      '/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.gz',
    #                      '/home/lab-heavy2/adna/bcf/shotgun/Yana_old/Yana_old_chr_notation_chr6.bcf')
    # run_only_concordance('/home/lab-heavy2/adna/output/simulated_SC_sorted/simulated_SC_non_MHC_region.bcf',
    #                      'chr6',
    #                      '/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.gz',
    #                      '/home/lab-heavy2/adna/bcf/shotgun/Yana_old/Yana_old_chr_notation_chr6.bcf')
    plot(data=os.path.join('/home/lab-heavy2/adna/output/simulated_SC_sorted/GLIMPSE_concordance_simulated_SC_MHC_region/output.rsquare.grp.txt.gz'),
         output=os.path.join('/home/lab-heavy2/adna/output/simulated_SC_sorted/accplot_MHC_region.png'))
    plot(data=os.path.join(
        '/home/lab-heavy2/adna/output/simulated_SC_sorted/GLIMPSE_concordance_simulated_SC_non_MHC_region/output.rsquare.grp.txt.gz'),
         output=os.path.join('/home/lab-heavy2/adna/output/simulated_SC_sorted/accplot_non_MHC_region.png'))
    # samples = ['Yana_old', 'VLASA32', 'PES001', 'I5241']
    # file_path_format = ('/home/lab-heavy2/adna/conclusions/heterozygosity/{sample_name}/'
    #                     'sample/{area}/{sample_name}_chr_notation_chr6{sites}.vcf.gz')
    # results = {}
    # for sample in samples:
    #     results[sample] = {}
    #     for area in [AREA_ALL, AREA_MHC, AREA_EXC_MHC]:
    #         results[sample][area] = {}
    #         for sites in ['all_sites', INC_1240K_NAME_ADDITION, EXC_1240K_NAME_ADDITION]:
    #             results[sample][area][sites] = {}
    #             counts = analyze_snp_depths_in_vcf(file_path_format.format(sample_name=sample,
    #                                                                        area=area,
    #                                                                        sites='' if sites == 'all_sites' else '_' + sites),
    #                                                4)
    #             results[sample][area][sites]['total'] = counts[0]
    #             results[sample][area][sites]['no_dp'] = counts[1]
    #             results[sample][area][sites]['dp_4_and_above'] = counts[2]
    #
    # with open('/home/lab-heavy2/adna/conclusions/sites_counts.json', 'w') as f:
    #     json.dump(results, f)

    # depths = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02]
    # downsample('/home/lab-heavy2/adna/bam/1240K/I5241/I5241_chr_notation_chr6.bam', depths, 0.265)
    # general_path = '/home/lab-heavy2/adna/output/Yana_old_chr_notation_chr6_subsampled_{level}x/renamed_ligated.bcf'
    # chr6_vcf_path = '/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/entire_chr6/Yana_old_chr_notation_chr6_subsampled_{level}x.vcf.gz'
    # MHC_vcf_path = '/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/MHC_region/Yana_old_chr_notation_chr6_subsampled_{level}x_MHC.vcf.gz'
    # exclude_MHC_vcf_path = '/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/chr6_excluding_MHC_region/Yana_old_chr_notation_chr6_subsampled_{level}x_exclude_MHC.vcf.gz'
    # for level in ['01', '025', '05', '075', '1', '2']:
    #     # create_vcf_from_bam(input_bam=general_path.format(level=level),
    #     #                     tsv='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.tsv.gz',
    #     #                     truth_fasta='/home/lab-heavy2/adna/reference/hg19.fa',
    #     #                     reference_sites='/home/lab-heavy2/adna/reference/hgdp1kgp/chr6/chr6.hg19.normalized.snps.sorted.chr6.fill_tags.sites.vcf.gz',
    #     #                     chromosome='chr6',
    #     #                     output_path=chr6_vcf_path.format(level=level))
    #     subprocess.run(['bcftools', 'view', general_path.format(level=level),
    #                     '-Oz', '-o', chr6_vcf_path.format(level=level)],
    #                    capture_output=True, text=True)
    #     subprocess.run(['tabix', '-p', 'vcf', chr6_vcf_path.format(level=level)])
    #     split_chr6_vcf(original_vcf_path=chr6_vcf_path.format(level=level),
    #                    mhc_region_file_path=MHC_vcf_path.format(level=level),
    #                    chr6_exclude_mhc_file_path=exclude_MHC_vcf_path.format(level=level))

    # filter_1240k_sites_from_all_vcf_files_in_dir('/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data')
