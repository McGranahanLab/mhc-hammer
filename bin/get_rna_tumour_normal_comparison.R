suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--allele',  nargs=1,
                    help='Allele',
                    required=TRUE)
parser$add_argument('--sample_name',  nargs=1,
                    help='Sample name',
                    required=TRUE)
parser$add_argument('--patient',  nargs=1,
                    help='Patient',
                    required=TRUE)

parser$add_argument('--allele_snp_bed',  nargs=1,
                    help='Allele1 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)

parser$add_argument('--allele_gl_coverage',  nargs=1,
                    help='Allele germline coverage file',
                    required=TRUE)
parser$add_argument('--allele_tumour_coverage',  nargs=1,
                    help='Allele tumour coverage file',
                    required=TRUE)

parser$add_argument('--tumour_library_size_path',  nargs=1,
                    help='Path to tumour flagstat csv',
                    required=TRUE)
parser$add_argument('--gl_library_size_path',  nargs=1,
                    help='Path to germline flagstat csv',
                    required=TRUE)

parser$add_argument('--output_path',  nargs=1,
                    help='Full path to output csv',
                    required=TRUE)
parser$add_argument('--plots_prefix',  nargs=1,
                    help='Full prefix to pdf plots',
                    required=TRUE)
parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

args <- parser$parse_args()

allele_name <- args$allele
sample_name <- args$sample_name
patient <- args$patient

allele_snp_bed <- args$allele_snp_bed

allele_gl_coverage_file <- args$allele_gl_coverage
allele_tumour_coverage_file <- args$allele_tumour_coverage

output_path <- args$output_path
plots_prefix <- args$plots_prefix

tumour_library_size_path <- args$tumour_library_size_path
gl_library_size_path <- args$gl_library_size_path

scripts_dir <- args$scripts_dir
source(paste0(scripts_dir, "mhc_hammer_functions.R"))

tumour_library_size <- as.numeric(readLines(tumour_library_size_path))
gl_library_size <- as.numeric(readLines(gl_library_size_path))

# perform comparison
tumour_normal_comparison <- calculate_tumour_normal_comparison(allele_name,
                                                               gl_library_size = gl_library_size,
                                                               tumour_library_size = tumour_library_size,
                                                               allele_mismatch_bed_file = allele_snp_bed,
                                                               allele_gl_coverage_file = allele_gl_coverage_file,
                                                               allele_tumour_coverage_file = allele_tumour_coverage_file,
                                                               make_plot = TRUE)

# save output
tumor_normal_output <- data.table(patient = patient,
                                  sample_name = sample_name,
                                  allele = allele_name, 
                                  paired_t_test = tumour_normal_comparison$paired_t_test,
                                  paired_wilcoxon_test = tumour_normal_comparison$paired_wilcoxon_test,
                                  median_tumour_dp = tumour_normal_comparison$median_tumour_rna_dp,
                                  median_normal_dp = tumour_normal_comparison$median_normal_rna_dp,
                                  normal_n_snps_with_coverage = tumour_normal_comparison$normal_n_snps_with_coverage,
                                  tumour_n_snps_with_coverage = tumour_normal_comparison$tumour_n_snps_with_coverage)

fwrite(tumor_normal_output, file = output_path)

# save figures

pdf(paste0(plots_prefix, "_boxplot.pdf"), width = 5)
print(tumour_normal_comparison$boxplot_plot)
dev.off()

pdf(paste0(plots_prefix, "_coverage.pdf"), width = 5, height = 2.5)
print(tumour_normal_comparison$coverage_plot)
dev.off()


