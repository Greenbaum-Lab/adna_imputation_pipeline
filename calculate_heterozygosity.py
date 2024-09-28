import pandas as pd
import numpy as np
import pysam
import os
import re
import argparse

REGIONS_POS = {'MHC': (28477797, 33448354),
               'MHC_CLASS_1': (29570015, 31478901),
               'MHC_CLASS_3': (31486754, 32374958),
               'MHC_CLASS_2': (32407666, 33377673)}
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
ALL_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'total_count', 'het_count']
TOT_BED_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'total_count']
HET_BED_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'het_count']
LOAD_COLUMNS = ['chromosome', 'start_pos', 'end_pos', 'count']
CSV_COLUMNS = ['sample_id', 'coverage', 'wg_het', 'MHC_class_1_het',
               'MHC_class_2_het', 'MHC_class_3_het', 'MHC_class_1_normalized_het',
               'MHC_class_2_normalized_het', 'MHC_class_3_normalized_het']
MIN_DP = 4
MIN_PROB = 0.9
OUTPUT_DIR = '/home/lab-heavy2/adna/output'
DEFAULT_COV_THRESHOLD = 1000
DEFAULT_SNP_THRESHOLD = 0
HLA_GENES_SITES_FILE = '/home/lab-heavy2/adna/hla_genes_sites.csv'


def calculate_heterozygosity_sites(input_file, is_imputed, sites_file,
                                   output_dir, chromosome, min_dp=MIN_DP,
                                   min_prob=MIN_PROB):
    """
    This function determines for each site in the sites_file whether in each
    sample it's homozygous, heterozygous, or has no data, and saves the results
    in two BEDGRAPH files.
    According to the is_imputed flag, there will be a different filter for sites
    that will be included in the calculation.

    :param input_file: To calculate heterozygosity level of
    :param is_imputed: whether the VCF file is an original sample or imputed data
    :param sites_file: File with relevant sites locations
    :param output_dir: dir to save output file at
    :param chromosome: to run calculations for
    :param min_dp: Minimum sequencing depth for original samples
    :param min_prob: Minimum probability for imputed data
    """
    bcf = pysam.VariantFile(input_file)
    sites = pd.read_csv(sites_file)
    sites = sites[['chromosome', 'position']]
    sites.rename(columns={'position': 'start_pos'}, inplace=True)
    sites['end_pos'] = sites['start_pos'] + 1
    sites['data'] = 0

    min_pos = sites.loc[0, 'start_pos']
    max_pos = sites.loc[len(sites) - 1, 'end_pos']

    current_line_index = 0

    for record in bcf.fetch(chromosome, min_pos, max_pos):
        while record.pos > sites.iloc[current_line_index]['start_pos']:
            # Mark -1 for no data for the site
            sites.loc[current_line_index, 'data'] = -1
            current_line_index += 1
        # If the site is not relevant based on sites_file we ignore it
        if record.pos != sites.iloc[current_line_index]['start_pos']:
            continue

        sample = record.samples[0]  # Access the first sample
        if sample.alleles and all(allele is not None for allele in sample.alleles):
            gt = sample['GT']
            dp = sample.get('DP')
            gp = sample.get('GP')
            if ((not is_imputed and gt and None not in gt and dp >= min_dp) or
                    (is_imputed and gt and None not in gt and max(gp) >= min_prob)):  # Check for valid genotype data
                is_het = len(set(gt)) > 1
                if is_het:
                    sites.loc[current_line_index, 'data'] = 1
                else:
                    sites.loc[current_line_index, 'data'] = 0
            else:
                sites.loc[current_line_index, 'data'] = -1
        else:
            sites.loc[current_line_index, 'data'] = -1

        current_line_index += 1

    while current_line_index < len(sites):
        sites.loc[current_line_index, 'data'] = -1
        current_line_index += 1

    bcf.close()

    output_file_name = os.path.basename(input_file).replace('.bcf',
                                                            '_sites_data.bg')
    sites.to_csv(os.path.join(output_dir, output_file_name),
                 sep='\t',
                 header=False,
                 index=False)


def calculate_heterozygosity_genes(input_file, is_imputed, genes_locations_file,
                                   output_dir, chromosome, maf_file, min_dp=MIN_DP,
                                   min_prob=MIN_PROB):
    """
    This function counts the total and heterozygous sites of a given BCF file
    in given genes positions, and saves the results in two BEDGRAPH files.
    According to the is_imputed flag, there will be a different filter for sites
    that will be included in the calculation.

    :param input_file: To calculate heterozygosity level of
    :param is_imputed: whether the VCF file is an original sample or imputed data
    :param genes_locations_file: file with gene names and positions
           (all in the same chromosome)
    :param output_dir: dir to save output file at
    :param chromosome: to run calculations for
    :param maf_file: list of positions in chromosome filtered based on MAF
    :param min_dp: Minimum sequencing depth for original samples
    :param min_prob: Minimum probability for imputed data
        """
    bcf = pysam.VariantFile(input_file)
    maf_include = pd.DataFrame(maf_file, delimiter=' ', names=['pos'])['pos'].values
    genes_data = pd.read_csv(genes_locations_file)
    genes_data.columns = ['chromosome', 'start_pos', 'end_pos']
    genes_data['chromosome'] = chromosome
    genes_data['total_count'] = 0
    genes_data['het_count'] = 0

    min_pos = genes_data.loc[0, 'start_pos']
    max_pos = genes_data.loc[len(genes_data) - 1, 'end_pos']

    for record in bcf.fetch(chromosome, min_pos, max_pos):
        sample = record.samples[0]  # Access the first sample
        if sample.alleles:  # Check if alleles information is available
            gt = sample['GT']
            dp = sample.get('DP')
            gp = sample.get('GP')
            if ((not is_imputed and gt and None not in gt and dp >= min_dp) or
                    (is_imputed and gt and None not in gt and max(gp) >= min_prob and record.pos in maf_include)):  # Check for valid genotype data
                is_het = len(set(gt)) > 1
                # Check if record overlaps any gene
                for idx, row in genes_data.iterrows():
                    if not (record.start >= row['end_pos'] or record.stop < row['start_pos']):
                        genes_data.loc[idx, 'total_count'] += 1
                        if is_het:
                            genes_data.loc[idx, 'het_count'] += 1

    bcf.close()

    tot_data = genes_data[TOT_BED_COLUMNS]
    tot_output_file_name = os.path.basename(input_file).replace('.bcf',
                                                                f'_total_counts_hla_genes.bg')
    tot_data.to_csv(os.path.join(output_dir, tot_output_file_name),
                    sep='\t',
                    header=False,
                    index=False)

    het_bins_data = genes_data[HET_BED_COLUMNS]
    het_output_file_name = os.path.basename(input_file).replace('.bcf',
                                                                f'_het_counts_hla_genes.bg')
    het_bins_data.to_csv(os.path.join(output_dir, het_output_file_name),
                         sep='\t',
                         header=False,
                         index=False)


def calculate_heterozygosity_bins(input_file, is_imputed, bin_size,
                                  output_dir, chromosome, maf_file, min_dp=MIN_DP,
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
    :param maf_file: list of positions in chromosome filtered based on MAF
    :param min_dp: Minimum sequencing depth for original samples
    :param min_prob: Minimum probability for imputed data
    """
    bcf = pysam.VariantFile(input_file)
    bins_data = pd.DataFrame(columns=ALL_COLUMNS)
    maf_include = pd.DataFrame(maf_file, delimiter=' ', names=['pos'])['pos'].values

    # Initiate bin data
    chrom_len = CHROM_LENGTHS[chromosome]
    start = 0
    while start < chrom_len:
        end = min(start + bin_size, chrom_len + 1)
        bins_data.loc[len(bins_data.index)] = [chromosome, start, end, 0, 0]
        start = end

    for record in bcf:
        sample = record.samples[0]  # Access the first sample
        if sample.alleles:  # Check if alleles information is available
            gt = sample['GT']
            dp = sample.get('DP')
            gp = sample.get('GP')
            if ((not is_imputed and gt and None not in gt and dp >= min_dp) or
                    (is_imputed and gt and None not in gt and max(gp) >= min_prob and record.pos in maf_include)): # Check for valid genotype data
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


def load_bed(bed_filename):
    data = pd.read_csv(bed_filename, delimiter='\t', names=LOAD_COLUMNS)
    return data


def calc_het_in_range(het_counts, total_counts, start_pos, end_pos):
    """
    This function calculates the heterozygosity in a given range.

    :param het_counts: pandas DataFrame that contains heterozygous sites counts
    :param total_counts: pandas DataFrame that contains total sites counts
    :param start_pos: Start position of the range.
    :param end_pos: End position of the range.
    :return: Heterozygosity in the given range.
    """
    in_range = np.all([total_counts['end_pos'].values > start_pos,
                       total_counts['start_pos'] < end_pos], axis=0)
    return het_counts[in_range]['count'].sum() / total_counts[in_range]['count'].sum()


def calc_het_outside_range(het_counts, total_counts, start_pos, end_pos):
    """
        This function calculates the heterozygosity outside the given range.

        :param het_counts: pandas DataFrame that contains heterozygous sites counts
        :param total_counts: pandas DataFrame that contains total sites counts
        :param start_pos: Start position of the range.
        :param end_pos: End position of the range.
        :return: Heterozygosity in the given range.
        """
    outside_range = np.any([total_counts['start_pos'].values > end_pos,
                            total_counts['end_pos'] < start_pos], axis=0)
    return het_counts[outside_range]['count'].sum() / total_counts[outside_range]['count'].sum()


def calculate_het_all_mhc_regions(het_counts, total_counts):
    """
    This method calculates the heterozygosity levels of4 MHC related regions:
        * MHC class 1
        * MHC class 2
        * MHC class 3
        * Non-MHC
    :param het_counts:
    :param total_counts:
    :return: tuple of length 4 with the heterozygosity levels of the regions
             in the order they are mentioned above.
    """
    return (calc_het_in_range(het_counts,
                              total_counts,
                              REGIONS_POS['MHC_CLASS_1'][0],
                              REGIONS_POS['MHC_CLASS_1'][1]),
            calc_het_in_range(het_counts,
                              total_counts,
                              REGIONS_POS['MHC_CLASS_2'][0],
                              REGIONS_POS['MHC_CLASS_2'][1]),
            calc_het_in_range(het_counts,
                              total_counts,
                              REGIONS_POS['MHC_CLASS_3'][0],
                              REGIONS_POS['MHC_CLASS_3'][1]),
            calc_het_outside_range(het_counts,
                                   total_counts,
                                   REGIONS_POS['MHC'][0],
                                   REGIONS_POS['MHC'][1]))


def find_call_coverage(total_counts, snp_threshold):
    """
    This method finds the sample's coverage (total snps count).

    :param total_counts: snps counts per bins
    :param snp_threshold:
    :return: total snp count
    """
    return np.sum(total_counts['count'].values > snp_threshold)


def create_het_csv(called_samples_pairs, imputed_samples_pairs, snp_threshold, coverage_threshold, output_dir):
    """
    This method created heterozygosity CSV for all given samples.
    If the coverage of a called sample is lower than the given threshold, the
    CSV records for it will contain np.nan (empty values)

    :param called_samples_pairs:
    :param imputed_samples_pairs:
    :param snp_threshold:
    :param coverage_threshold: to filter samples
    :param output_dir: path to save output CSV at
    """
    call_csv_data = []
    phase_csv_data = []

    for sample_name in called_samples_pairs:
        called_het_counts_file = called_samples_pairs[sample_name]['het_counts']
        called_total_counts_file = called_samples_pairs[sample_name]['total_counts']
        imputed_het_counts_file = imputed_samples_pairs[sample_name]['het_counts']
        imputed_total_counts_file = imputed_samples_pairs[sample_name]['total_counts']

        if (not os.path.isfile(called_het_counts_file) or
                not os.path.isfile(called_total_counts_file) or
                not os.path.isfile(imputed_het_counts_file) or
                not os.path.isfile(imputed_total_counts_file)):
            raise FileNotFoundError('No heterozygosity counts files found for the given sample')

        called_het_counts = load_bed(called_het_counts_file)
        called_total_counts = load_bed(called_total_counts_file)
        imputed_het_counts = load_bed(imputed_het_counts_file)
        imputed_total_counts = load_bed(imputed_total_counts_file)

        # Check if the sample should be filtered
        coverage = find_call_coverage(called_total_counts, snp_threshold)
        if coverage < coverage_threshold:
            # call_csv_data.append([sample_name, coverage] + [np.nan] * 7)
            # phase_csv_data.append([sample_name, coverage] + [np.nan] * 7)
            continue

        called_class_1_het, called_class_2_het, called_class_3_het, called_non_mhc_het = (
            calculate_het_all_mhc_regions(called_het_counts, called_total_counts))
        call_csv_data.append([sample_name, coverage, called_non_mhc_het,
                              called_class_1_het, called_class_2_het,
                              called_class_3_het,
                              called_class_1_het - called_non_mhc_het,
                              called_class_2_het - called_non_mhc_het,
                              called_class_3_het - called_non_mhc_het])

        imputed_class_1_het, imputed_class_2_het, imputed_class_3_het, imputed_non_mhc_het = (
            calculate_het_all_mhc_regions(imputed_het_counts, imputed_total_counts))
        phase_csv_data.append([sample_name, coverage, imputed_non_mhc_het,
                               imputed_class_1_het, imputed_class_2_het,
                               imputed_class_3_het,
                               imputed_class_1_het - imputed_non_mhc_het,
                               imputed_class_2_het - imputed_non_mhc_het,
                               imputed_class_3_het - imputed_non_mhc_het])

    called_df = pd.DataFrame(call_csv_data)
    called_df.columns = CSV_COLUMNS
    called_df.to_csv(os.path.join(output_dir, 'called_samples_het.csv'),
                     header=True, index=False)

    imputed_df = pd.DataFrame(phase_csv_data)
    imputed_df.columns = CSV_COLUMNS
    imputed_df.to_csv(os.path.join(output_dir, 'imputed_samples_het.csv'),
                      header=True, index=False)


def find_file_pairs(input_dir):
    # Regular expression to parse the file names
    pattern = re.compile(r'^(.*)_((?:het_counts)|(?:total_counts))_bin_size_\d+\.bg$')

    # Dictionary to hold the file pairs
    file_pairs = {}

    # List all files in the directory
    for filename in os.listdir(input_dir):
        # Match the file name to the pattern
        match = pattern.match(filename)
        if match:
            sample_name, count_type = match.groups()[0:2]
            if sample_name not in file_pairs:
                file_pairs[sample_name] = {}
            if count_type == 'het_counts':
                file_pairs[sample_name]['het_counts'] = os.path.join(input_dir, filename)
            elif count_type == 'total_counts':
                file_pairs[sample_name]['total_counts'] = os.path.join(input_dir, filename)

    # Filter out any incomplete pairs
    complete_pairs = {k: v for k, v in file_pairs.items() if 'het_counts' in v and 'total_counts' in v}

    return complete_pairs


def create_csv_for_all_samples_in_dir(input_dir,
                                      coverage_threshold=DEFAULT_COV_THRESHOLD,
                                      snp_threshold=DEFAULT_SNP_THRESHOLD):
    """

    :param input_dir:
    :param coverage_threshold: minimal number of bins with data in sample
    :param snp_threshold: minimum snps in bin for it to be counted
    :return:
    """
    called_samples = find_file_pairs(os.path.join(input_dir, 'called'))
    imputed_samples = find_file_pairs(os.path.join(input_dir, 'imputed'))

    create_het_csv(called_samples, imputed_samples, snp_threshold,
                   coverage_threshold, input_dir)
