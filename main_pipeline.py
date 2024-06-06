# Imports
import os.path
import subprocess
import argparse

from glimpse_wrapper import GlimpseWrapper
from utils import create_bcf_from_bam
from calculate_heterozygosity import (calculate_heterozygosity_bins,
                                      calculate_heterozygosity_genes)

# Consts
GLIMPSE_DIR = '/home/lab-heavy2/Documents/yuvalt/glimpse'
WD = '/home/lab-heavy2/adna'
TRUTH_FASTA = '/home/lab-heavy2/adna/reference/hg19.fa'
DEFAULT_BIN_SIZE = 10000
DEFAULT_GENES_POS_FILE = '~/keith/adna_db/hla_genes_positions'


def call(input_bam, output_path, chromosome, tsv, reference_sites,
         truth_fasta=TRUTH_FASTA):
    """
    This method converts the given BAM file into a BCF file located and named
    in the given output path.

    :param input_bam: BAM to convert
    :param output_path: path to output file
    :param chromosome:
    :param tsv:
    :param reference_sites:
    :param truth_fasta:
    """
    if not reference_sites or not output_path:
        raise TypeError('call mode requires the following arguments: '
                        'input_bam, output_path, chromosome, reference_sites'
                        ' (, tsv)')
    create_bcf_from_bam(input_bam=input_bam,
                        chromosome=chromosome,
                        tsv=tsv,
                        truth_fasta=truth_fasta,
                        reference_sites=reference_sites,
                        output_path=output_path)


def phase(input_bam, reference, chromosome, reference_sites=None,
          output_path=None, glimpse_dir=GLIMPSE_DIR, wd=WD,
          truth_fasta=TRUTH_FASTA):
    """
    This method runs imputation on the input BAM file.

    :param input_bam:
    :param output_path:
    :param reference:
    :param chromosome:
    :param reference_sites:
    :param tsv:
    :param glimpse_dir:
    :param wd:
    :param truth_fasta:
    :return:
    """
    if not reference_sites:
        raise TypeError('phase mode requires the following arguments: '
                        'input_bam, reference, chromosome (, reference_sites'
                        ', output_path)')
    g = GlimpseWrapper(glimpse_dir=glimpse_dir,
                       wd=wd,
                       input_bam=input_bam,
                       chromosome=chromosome,
                       reference=reference,
                       reference_sites=reference_sites,
                       truth_fasta=truth_fasta)
    g.run()

    if output_path:
        subprocess.run(['mv', g.get_output_path(), output_path],
                       capture_output=True, text=True)
        subprocess.run(['mv', g.get_output_path() + '.csi', output_path + '.csi'],
                       capture_output=True, text=True)

    g.delete_temp_files()


def het(input_bcf, output_dir, is_imputed, chromosome, bin_size=DEFAULT_BIN_SIZE,
        genes_pos_file=DEFAULT_GENES_POS_FILE):
    """
    This method runs heterozygosity sites counts on the given BCF file.

    :param input_bcf:
    :param output_dir:
    :param is_imputed:
    :param chromosome:
    :param bin_size:
    :param genes_pos_file:
    """
    calculate_heterozygosity_bins(input_bcf, is_imputed, bin_size, output_dir, chromosome)
    calculate_heterozygosity_genes(input_bcf, is_imputed, genes_pos_file, output_dir, chromosome)


def main():
    parser = argparse.ArgumentParser(description='Analyze input BAM file in different modes')
    parser.add_argument('mode', type=str, choices=['call', 'phase', 'het'],
                        help='Mode for running the pipeline')
    parser.add_argument('input_path', type=str,
                        help='Path to input file (BAM for "call" and "phase" modes, BCF for "het" mode')
    parser.add_argument('--output_path', type=str, help='Path to output file')
    parser.add_argument('--chromosome', type=str, required=True)
    parser.add_argument('--reference', type=str, help='Path to reference file')
    parser.add_argument('--tsv', type=str,
                        help='Path to tsv (used in formats conversion).'
                             ' Not required - can be created automatically'
                             ' given reference_sites file.')
    parser.add_argument('--reference_sites', type=str,
                        help='Path to reference sites file. Not required - '
                             'can be created automatically given a reference file.')
    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('--imputed', dest='is_imputed', action='store_true')
    feature_parser.add_argument('--non-imputed', dest='is_imputed', action='store_false')
    parser.set_defaults(is_imputed=None)
    parser.add_argument('--glimpse_dir', type=str,
                        help='Directory of GLIMPSE source code')
    parser.add_argument('--wd', type=str,
                        help='Directory to create GLIMPSE temp and output files at.')
    parser.add_argument('--truth_fasta', type=str,
                        help='The faidx-indexed reference file in the FASTA format.')
    parser.add_argument('--genes_pos_file', type=str,
                        help='Path to file containing the genes positions to '
                             'calculate heterozygosity for.')

    args = parser.parse_args()

    # Validate arguments groups
    if args.mode == 'call':
        if args.output_path is None or args.reference_sites is None:
            raise ValueError("call mode requires arguments: input_path, "
                             "output_path, chromosome, reference_sites, truth_fasta [,tsv]")
        call(args.input_path, args.output_path, args.chromosome,
             args.tsv, args.reference_sites, args.truth_fasta)
    elif args.mode == 'phase':
        if args.reference is None:
            raise ValueError("phase mode requires arguments: input_path, "
                             "chromosome, reference [, output_path, reference_sites,"
                             "glimpse_dir, wd, truth_fasta]")
        phase(args.input_path, args.reference, args.chromosome,
              args.reference_sites, args.output_path, args.glimpse_dir,
              args.wd, args.truth_fasta)
    else:
        if args.is_imputed is None:
            raise ValueError("het mode requires arguments: input_path, chromosome, imputed / non-imputed")
        if not args.output_path:
            het(args.input_path, os.path.dirname(args.input_path), args.is_imputed,
                args.chromosome, DEFAULT_BIN_SIZE, args.genes_pos_file)
        else:
            het(args.input_path, args.output_path, args.is_imputed,
                args.chromosome, DEFAULT_BIN_SIZE, args.genes_pos_file)


if __name__ == '__main__':
    main()
