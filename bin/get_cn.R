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

parser$add_argument('--allele1_gl_coverage_file',  nargs=1,
                    help='Allele 1 germline coverage file',
                    required=TRUE)
parser$add_argument('--allele2_gl_coverage_file',  nargs=1,
                    help='Allele 2 germline coverage file',
                    required=TRUE)

parser$add_argument('--allele1_tumour_coverage_file',  nargs=1,
                    help='Allele 1 tumour coverage file',
                    required=TRUE)
parser$add_argument('--allele2_tumour_coverage_file',  nargs=1,
                    help='Allele 2 tumour coverage file',
                    required=TRUE)

parser$add_argument('--allele1_snp_bed',  nargs=1,
                    help='Allele1 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)
parser$add_argument('--allele2_snp_bed',  nargs=1,
                    help='Allele2 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)

parser$add_argument('--purity',  nargs=1,
                    help='Purity (range 0-1)',
                    required=TRUE)
parser$add_argument('--ploidy',  nargs=1,
                    help='Ploidy',
                    required=TRUE)
parser$add_argument('--gtf_path',  nargs=1,
                    help='Path to GTF file',
                    required=TRUE)
parser$add_argument('--cn_output_path',  nargs=1,
                    help='Full path to output csv',
                    required=TRUE)
parser$add_argument('--cn_plots_prefix',  nargs=1,
                    help='Full prefix to pdf plots',
                    required=TRUE)

parser$add_argument('--tumour_library_size_path',  nargs=1,
                    help='Path to tumour library size',
                    required=TRUE)
parser$add_argument('--gl_library_size_path',  nargs=1,
                    help='Path to germline library size',
                    required=TRUE)

parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)


args <- parser$parse_args()

allele1 <- args$allele1
allele2 <- args$allele2

allele1_gl_coverage_file <- args$allele1_gl_coverage_file
allele2_gl_coverage_file <- args$allele2_gl_coverage_file

allele1_tumour_coverage_file <- args$allele1_tumour_coverage_file
allele2_tumour_coverage_file <- args$allele2_tumour_coverage_file

allele1_snp_bed <- args$allele1_snp_bed
allele2_snp_bed <- args$allele2_snp_bed

purity <- as.numeric(args$purity)
ploidy <- as.numeric(args$ploidy)
gtf_path <- args$gtf_path
cn_output_path <- args$cn_output_path
cn_plots_prefix <- args$cn_plots_prefix

tumour_library_size_path <- args$tumour_library_size_path
gl_library_size_path <- args$gl_library_size_path

scripts_dir <- args$scripts_dir

source(paste0(scripts_dir, "/mhc_hammer_functions.R"))

tumour_library_size <- as.numeric(readLines(tumour_library_size_path))
gl_library_size <- as.numeric(readLines(gl_library_size_path))

if(is.na(tumour_library_size) | is.na(gl_library_size)){
  stop("Either tumour or germline library size is NA")
}

cn <- calculate_cn(bin_size =  150, gamma = 1,
                   allele1_name = allele1,
                   allele2_name = allele2,
                   allele1_gl_coverage_file = allele1_gl_coverage_file,
                   allele1_tumour_coverage_file = allele1_tumour_coverage_file,
                   allele2_gl_coverage_file = allele2_gl_coverage_file,
                   allele2_tumour_coverage_file = allele2_tumour_coverage_file,
                   allele1_snp_bed_file = allele1_snp_bed,
                   allele2_snp_bed_file = allele2_snp_bed,
                   purity = purity,
                   ploidy = ploidy,
                   gl_library_size = gl_library_size,
                   tumour_library_size = tumour_library_size,
                   gtf_path = gtf_path,
                   make_plot = TRUE)

# save output
cn_output <- data.table(allele1 = allele1,
                        allele2 = allele2,
                        cn1 = cn$cn1,
                        cn1_lower = cn$cn1_lower,
                        cn1_upper = cn$cn1_upper,
                        cn2 = cn$cn2,
                        cn2_lower = cn$cn2_lower,
                        cn2_upper = cn$cn2_upper,
                        cn1_binned = cn$cn1_binned,
                        cn1_binned_lower = cn$cn1_binned_lower,
                        cn1_binned_upper = cn$cn1_binned_upper,
                        cn2_binned = cn$cn2_binned,
                        cn2_binned_lower = cn$cn2_binned_lower,
                        cn2_binned_upper = cn$cn2_binned_upper,
                        n_bins = cn$n_bins,
                        n_snps = cn$n_snps)

fwrite(cn_output, file = cn_output_path)

# save figures 

pdf(paste0(cn_plots_prefix, "_boxplot.pdf"), width = 12)
print(ggarrange(cn$cn_boxplot_plot, cn$cn_bin_boxplot_plot, ncol = 2))
dev.off()

pdf(paste0(cn_plots_prefix, ".pdf"), width = 12)
print(ggarrange(cn$allele1_cn_plot, cn$allele2_cn_plot, ncol = 1))
dev.off()

pdf(paste0(cn_plots_prefix, "_logr.pdf"), width = 12)
print(ggarrange(cn$allele1_logr_plot, cn$allele2_logr_plot, ncol = 1))
dev.off()

pdf(paste0(cn_plots_prefix, "_coverage.pdf"), width = 12)
print(ggarrange(cn$allele1_coverage_plot, cn$allele2_coverage_plot, ncol = 1))
dev.off()
