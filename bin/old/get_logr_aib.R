suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))

parser <- ArgumentParser()

parser$add_argument('--allele1',  nargs=1,
                    help='Allele 1',
                    required=TRUE)
parser$add_argument('--allele2',  nargs=1,
                    help='Allele2',
                    required=TRUE)

parser$add_argument('--allele1_gl_reads_count_once_coverage',  nargs=1,
                    help='Allele 1 germline coverage file',
                    required=TRUE)
parser$add_argument('--allele2_gl_reads_count_once_coverage',  nargs=1,
                    help='Allele 2 germline coverage file',
                    required=TRUE)

parser$add_argument('--allele1_tumour_reads_count_once_coverage',  nargs=1,
                    help='Allele 1 tumour coverage file',
                    required=TRUE)
parser$add_argument('--allele2_tumour_reads_count_once_coverage',  nargs=1,
                    help='Allele 2 tumour coverage file',
                    required=TRUE)

parser$add_argument('--allele1_snp_bed',  nargs=1,
                    help='Allele1 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)
parser$add_argument('--allele2_snp_bed',  nargs=1,
                    help='Allele2 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)

parser$add_argument('--tumour_library_size_path',  nargs=1,
                    help='Path to tumour library size',
                    required=TRUE)
parser$add_argument('--gl_library_size_path',  nargs=1,
                    help='Path to germline library size',
                    required=TRUE)

parser$add_argument('--logr_aib_output_path',  nargs=1,
                    help='Full path to output csv',
                    required=TRUE)
parser$add_argument('--logr_aib_plots_prefix',  nargs=1,
                    help='Full prefix to pdf plots',
                    required=TRUE)
parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)


args <- parser$parse_args()

allele1 <- args$allele1
allele2 <- args$allele2

allele1_gl_reads_count_once_coverage <- args$allele1_gl_reads_count_once_coverage
allele2_gl_reads_count_once_coverage <- args$allele2_gl_reads_count_once_coverage

allele1_tumour_reads_count_once_coverage <- args$allele1_tumour_reads_count_once_coverage
allele2_tumour_reads_count_once_coverage <- args$allele2_tumour_reads_count_once_coverage

allele1_snp_bed <- args$allele1_snp_bed
allele2_snp_bed <- args$allele2_snp_bed

tumour_library_size_path <- args$tumour_library_size_path
gl_library_size_path <- args$gl_library_size_path

logr_aib_output_path <- args$logr_aib_output_path
logr_aib_plots_prefix <- args$logr_aib_plots_prefix

scripts_dir <- args$scripts_dir

source(paste0(scripts_dir, "/mhc_hammer_functions.R"))

tumour_library_size <- as.numeric(readLines(tumour_library_size_path))
gl_library_size <- as.numeric(readLines(gl_library_size_path))

aib <- calculate_logr_aib(gl_library_size = gl_library_size,
                          tumour_library_size = tumour_library_size,
                          allele1_name = allele1,
                          allele2_name = allele2,
                          allele1_snp_bed_file = allele1_snp_bed,
                          allele2_snp_bed_file = allele2_snp_bed,
                          allele1_gl_reads_count_once_coverage = allele1_gl_reads_count_once_coverage,
                          allele2_gl_reads_count_once_coverage = allele2_gl_reads_count_once_coverage,
                          allele1_tumour_reads_count_once_coverage = allele1_tumour_reads_count_once_coverage,
                          allele2_tumour_reads_count_once_coverage = allele2_tumour_reads_count_once_coverage,
                          make_plot = TRUE)

# save output
aib_output <- data.table(allele1 = allele1,
                         allele2 = allele2,
                         n_snps = aib$n_snps,
                         paired_t_test = aib$paired_t_test,
                         paired_wilcoxon_test = aib$paired_wilcoxon_test)

fwrite(aib_output, file = logr_aib_output_path)

# save figures 

pdf(paste0(logr_aib_plots_prefix, "_boxplot.pdf"), width = 5)
print(aib$boxplot)
dev.off()

