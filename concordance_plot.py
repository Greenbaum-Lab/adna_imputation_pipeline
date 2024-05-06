# Imports
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import gzip


def plot(data='./GLIMPSE_concordance/output.rsquare.grp.txt.gz',
         output='accplot.png'):
    x_maf = []
    y_aggr = []

    plt.rc('axes', labelsize=20)

    with gzip.open(data, 'rt') as csvfile:
        plots = csv.reader(csvfile, delimiter=' ')
        for row in plots:
            y_aggr.append(float(row[4]))
            x_maf.append(float(row[2]))

    fig, axs = plt.subplots(1)

    axs.semilogx(x_maf, y_aggr, '--',  marker='o', lw=1, label='Aggregate r2',
                 markersize=10, alpha=0.8)

    mybins = [0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15,
              0.20, 0.3, 0.4, 0.50, 0.60]
    labels = ['0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5', '10', '15', '20',
              '30', '40', '50', '60']
    axs.set_xlabel('Minor allele frequency (%)')
    axs.set_ylabel('$r^2$ imputed vs true genotypes')
    axs.grid(linestyle='--')
    axs.legend(loc="lower right", prop={'size': 22})
    axs.set_ylim([0, 1.05])
    axs.set_xlim([0.0005, 0.6])
    axs.set_xticks(mybins)
    axs.set_xticklabels(labels)
    axs.minorticks_off()
    axs.minorticks_off()

    axs.set_title('Imputation accuracy', fontsize=20)
    fig.set_size_inches(16, 12)
    fig.savefig(output, dpi=300)


if __name__ == '__main__':
    for level in ['025', '05', '075', '1', '2']:
        plot(f'/home/lab-heavy2/adna/output/Yana_old_chr_notation_chr6_subsampled_{level}x/GLIMPSE_concordance_MHC_region/output.rsquare.grp.txt.gz',
             f'/home/lab-heavy2/adna/output/Yana_old_chr_notation_chr6_subsampled_{level}x/accplot_MHC_region.png')
        plot(f'/home/lab-heavy2/adna/output/Yana_old_chr_notation_chr6_subsampled_{level}x/GLIMPSE_concordance_chr6_excluding_MHC_region/output.rsquare.grp.txt.gz',
             f'/home/lab-heavy2/adna/output/Yana_old_chr_notation_chr6_subsampled_{level}x/accplot_chr6_excluding_MHC_region.png')
