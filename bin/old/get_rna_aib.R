suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--patient',  nargs=1,
                    help='Patient name',
                    required=TRUE)
parser$add_argument('--sample_name',  nargs=1,
                    help='Sample name',
                    required=TRUE)

parser$add_argument('--allele1',  nargs=1,
                    help='Allele 1',
                    required=TRUE)
parser$add_argument('--allele2',  nargs=1,
                    help='Allele2',
                    required=TRUE)

parser$add_argument('--allele1_reads_count_once_coverage',  nargs=1,
                    help='Allele 1 germline coverage file',
                    required=TRUE)
parser$add_argument('--allele2_reads_count_once_coverage',  nargs=1,
                    help='Allele 2 germline coverage file',
                    required=TRUE)

parser$add_argument('--allele1_snp_bed',  nargs=1,
                    help='Allele1 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)
parser$add_argument('--allele2_snp_bed',  nargs=1,
                    help='Allele2 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)

parser$add_argument('--aib_output_path',  nargs=1,
                    help='Full path to output csv',
                    required=TRUE)
parser$add_argument('--aib_plots_prefix',  nargs=1,
                    help='Full prefix to pdf plots',
                    required=TRUE)

parser$add_argument('--plot_title',  nargs=1,
                    help='Box plot title',
                    required=TRUE)

parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

args <- parser$parse_args()

patient <- args$patient
sample_name <- args$sample_name

allele1_name <- args$allele1
allele2_name <- args$allele2

allele1_mismatch_bed_file <- args$allele1_snp_bed
allele2_mismatch_bed_file <- args$allele2_snp_bed

allele1_unique_coverage_file <- args$allele1_reads_count_once_coverage
allele2_unique_coverage_file <- args$allele2_reads_count_once_coverage

aib_output_path <- args$aib_output_path
aib_plots_prefix <- args$aib_plots_prefix

plot_title <- args$plot_title

scripts_dir <- args$scripts_dir
source(paste0(scripts_dir, "mhc_hammer_functions.R"))

aib <- calculate_depth_aib(allele1_name = allele1_name,
                           allele2_name = allele2_name,
                           allele1_mismatch_bed_file = allele1_mismatch_bed_file,
                           allele2_mismatch_bed_file = allele2_mismatch_bed_file,
                           allele1_unique_coverage_file = allele1_unique_coverage_file,
                           allele2_unique_coverage_file = allele2_unique_coverage_file,
                           plot_title = plot_title,
                           make_plot = TRUE)

# save output
aib_output <- data.table(patient = patient,
                         sample_name = sample_name,
                         allele1 = allele1_name, 
                         allele2 = allele2_name, 
                         paired_t_test = aib$paired_t_test,
                         paired_wilcoxon_test = aib$paired_wilcoxon_test,
                         allele1_n_snps_with_coverage = aib$allele1_n_snps_with_coverage,
                         allele2_n_snps_with_coverage = aib$allele2_n_snps_with_coverage)

fwrite(aib_output, file = aib_output_path)

# save figures

pdf(paste0(aib_plots_prefix, "_boxplot.pdf"), width = 5)
print(aib$boxplot)
dev.off()
