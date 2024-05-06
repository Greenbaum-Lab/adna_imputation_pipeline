import pandas as pd
import pysam
import os
import json


MHC_START_POS = 28477797
MHC_END_POS = 33448354
CHROM_LENGTHS = {'chr1': 249250621,
                 'chr2': 243199373,
                 'chr3': 198022430,
                 'chr4': 191154276,
                 'chr5': 180915260,
                 'chr6': 171115067,
                 'chr7': 159138663,
                 'chr8': 146364022,
                 'chr9': 141213431,
                 'chr10': 135534747,
                 'chr11': 135006516,
                 'chr12': 133851895,
                 'chr13': 115169878,
                 'chr14': 107349540,
                 'chr15': 102531392,
                 'chr16': 90354753,
                 'chr17': 81195210,
                 'chr18': 78077248,
                 'chr19': 59128983,
                 'chr20': 63025520,
                 'chr21': 48129895,
                 'chr22': 51304566,
                 'chrX': 155270560,
                 'chrY': 59373566}
COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'total_count', 'het_count']
TOT_BED_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'total_count']
HET_BED_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'het_count']
MIN_DP = 4
MIN_PROB = 0.9


def calculate_heterozygosity_bins(input_file, is_imputed, bin_size,
                                  output_dir, chromosome,  min_dp=MIN_DP,
                                  min_prob=MIN_PROB):
    """
    This function counts the total and heterozygous sites of a given BCF file
    in division to bins of bin_size, and saves the results in two BEDGRAPH
    files.
    According to the is_imputed flag, there will be a different filter for sites
    that will be included in the calculation.

    :param input_file: To calculate heterozygosity level of
    :param is_imputed: whether the VCF file is an original sample or imputed data
    :param bin_size: size of bins to divide the calculation for
    :param output_dir: dir to save output file at
    :param chromosome: to run calculations for
    :param min_dp: Minimum sequencing depth for original samples
    :param min_prob: Minimum probability for imputed data
    """
    bcf = pysam.VariantFile(input_file)
    bins_data = pd.DataFrame(columns=COLUMNS)

    # Initiate bin data
    chrom_len = CHROM_LENGTHS[chromosome]
    start = 0
    while start < chrom_len:
        end = min(start + bin_size + 1, chrom_len + 1)
        bins_data.loc[len(bins_data.index)] = [chromosome, start, end, 0, 0]
        start = end

    for record in bcf:
        sample = record.samples[0]  # Access the first sample
        if sample.alleles:  # Check if alleles information is available
            gt = sample['GT']
            dp = sample.get('DP')
            gp = sample.get('GP')
            if ((not is_imputed and gt and None not in gt and dp >= min_dp) or
                    (is_imputed and gt and None not in gt and max(gp) >= min_prob)):  # Check for valid genotype data
                bin_index = record.pos // bin_size
                bins_data.loc[bin_index, 'total_count'] += 1
                if len(set(gt)) > 1:  # Heterozygous if more than one allele in the genotype
                    bins_data.loc[bin_index, 'het_count'] += 1

    bcf.close()

    tot_bins_data = bins_data[TOT_BED_COLUMNS]
    tot_output_file_name = os.path.basename(input_file).replace('.bcf',
                                                                f'_total_counts_bin_size_{bin_size}.bg')
    tot_bins_data.to_csv(os.path.join(output_dir, tot_output_file_name),
                         sep='\t',
                         header=False,
                         index=False)

    het_bins_data = bins_data[HET_BED_COLUMNS]
    het_output_file_name = os.path.basename(input_file).replace('.bcf',
                                                                f'_het_counts_bin_size_{bin_size}.bg')
    het_bins_data.to_csv(os.path.join(output_dir, het_output_file_name),
                         sep='\t',
                         header=False,
                         index=False)


def calculate_heterozygosity_whole_file(input_file, is_imputed, min_dp=MIN_DP, min_prob=MIN_PROB):
    """
    This function calculates a VCF file's heterozygosity level.
    According to the is_imputed flag, there will be a different filter for sites
    that will be included in the calculation.

    :param input_file: To calculate heterozygosity level of
    :param is_imputed: whether the VCF file is an original sample or imputed data
    :param min_dp: Minimum sequencing depth for original samples
    :param min_prob: Minimum probability for imputed data
    :return: Heterozygosity level, heterozygous sites count, total sites count
    """
    vcf = pysam.VariantFile(input_file)
    heterozygous_count = 0
    total_count = 0

    for record in vcf:
        sample = record.samples[0]  # Access the first sample
        if sample.alleles:  # Check if alleles information is available
            gt = sample['GT']
            dp = sample.get('DP')
            gp = sample.get('GP')
            if ((not is_imputed and gt and None not in gt and dp >= min_dp) or
                    (is_imputed and gt and None not in gt and max(gp) >= min_prob)):  # Check for valid genotype data
                total_count += 1
                if len(set(gt)) > 1:  # Heterozygous if more than one allele in the genotype
                    heterozygous_count += 1

    vcf.close()

    # Calculate and return heterozygosity
    if total_count > 0:
        return heterozygous_count / total_count, heterozygous_count, total_count


if __name__ == '__main__':
    print(calculate_heterozygosity_whole_file('/home/lab-heavy2/adna/conclusions/heterozygosity/Yana_old/imputed_data/entire_chr6/Yana_old_chr_notation_chr6.vcf.gz',
                                              is_imputed=True, min_prob=0.9))
    # samples = {'Yana_old': '/home/lab-heavy2/adna/bam/shotgun/Yana_old/Yana_old_chr_notation_chr6.bam',
    #            'VLASA32': '/home/lab-heavy2/adna/bam/shotgun/VLASA32/VLASA32_chr_notation_chr6.bam',
    #            'PES001': '/home/lab-heavy2/adna/bam/shotgun/PES001/PES001_chr_notation_chr6.bam',
    #            'I5241': '/home/lab-heavy2/adna/bam/1240K/I5241/I5241_chr_notation_chr6.bam'}
    # levels = [0.1, 0.25, 0.5, 0.75, 1, 2]
    # results = {'entire_chr6': {},
    #            'MHC_region': {},
    #            'chr6_excluding_MHC_region': {}}
    # chr6_vcf_path = ('/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/'
    #                  'entire_chr6/Yana_old_chr_notation_chr6_subsampled_{level}x_exclude_1240K_sites.vcf.gz')
    # MHC_vcf_path = ('/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/'
    #                 'MHC_region/Yana_old_chr_notation_chr6_subsampled_{level}x_MHC_exclude_1240K_sites.vcf.gz')
    # exclude_MHC_vcf_path = ('/home/lab-heavy2/adna/conclusions/heterozygosity/imputed_data/'
    #                         'chr6_excluding_MHC_region/'
    #                         'Yana_old_chr_notation_chr6_subsampled_{level}x_exclude_MHC_exclude_1240K_sites.vcf.gz')
    # for level in levels:
    #     results['entire_chr6'][level] = (
    #         calculate_heterozygosity(chr6_vcf_path.format(level=str(level).replace('.', ''))))
    #     results['MHC_region'][level] = (
    #         calculate_heterozygosity(MHC_vcf_path.format(level=str(level).replace('.', ''))))
    #     results['chr6_excluding_MHC_region'][level] = (
    #         calculate_heterozygosity(exclude_MHC_vcf_path.format(level=str(level).replace('.', ''))))
    #
    # with open('/home/lab-heavy2/adna/conclusions/Yana_old_exclude_1240K_sites_imputed_data_downsampling_levels_to_heterozygosity.json', 'w') as f:
    #     json.dump(results, f)

