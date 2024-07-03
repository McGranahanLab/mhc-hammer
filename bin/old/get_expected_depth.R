suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

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

parser$add_argument('--allele1_snp_bed',  nargs=1,
                    help='Allele1 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)
parser$add_argument('--allele2_snp_bed',  nargs=1,
                    help='Allele2 SNP bed file, only including SNPs that pass minimum coverage in both alleles',
                    required=TRUE)

parser$add_argument('--purity',  nargs=1,
                    help='Purity (range 0-1)',
                    required=TRUE)

parser$add_argument('--expected_depth_output_path',  nargs=1,
                    help='Full path to output csv',
                    required=TRUE)
parser$add_argument('--expected_depth_plots_prefix',  nargs=1,
                    help='Full prefix to pdf plots',
                    required=TRUE)

parser$add_argument('--tumour_library_size_path',  nargs=1,
                    help='Path to tumour flagstat csv',
                    required=TRUE)
parser$add_argument('--gl_library_size_path',  nargs=1,
                    help='Path to germline flagstat csv',
                    required=TRUE)

parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

args <- parser$parse_args()

allele1 <- args$allele1
allele2 <- args$allele2

allele1_gl_coverage_file <- args$allele1_gl_coverage_file
allele2_gl_coverage_file <- args$allele2_gl_coverage_file

allele1_snp_bed <- args$allele1_snp_bed
allele2_snp_bed <- args$allele2_snp_bed
gtf_path <- args$gtf_path

purity <- as.numeric(args$purity)
expected_depth_output_path <- args$expected_depth_output_path
expected_depth_plots_prefix <- args$expected_depth_plots_prefix

tumour_library_size_path <- args$tumour_library_size_path
gl_library_size_path <- args$gl_library_size_path

scripts_dir <- args$scripts_dir
source(paste0(scripts_dir, "mhc_hammer_functions.R"))

tumour_library_size <- as.numeric(readLines(tumour_library_size_path))
gl_library_size <- as.numeric(readLines(gl_library_size_path))

if(is.na(tumour_library_size) | is.na(gl_library_size)){
  stop("Either tumour or germline library size is NA")
}

expected_depth_all_snps <- calculate_exp_dp(allele1_name = allele1,
                                   allele2_name = allele2,
                                   allele1_gl_coverage_file = allele1_gl_coverage_file,
                                   allele2_gl_coverage_file = allele2_gl_coverage_file,
                                   allele1_snp_bed_file = allele1_snp_bed,
                                   allele2_snp_bed_file = allele2_snp_bed,
                                   purity = purity,
                                   gl_library_size = gl_library_size,
                                   tumour_library_size = tumour_library_size)

# save output
expected_depth_output <- data.table(allele1 = allele1,
                                    allele2 = allele2,
                                    allele1_exp_dp = expected_depth_all_snps$allele1_exp_dp,
                                    allele2_exp_dp = expected_depth_all_snps$allele2_exp_dp)

fwrite(expected_depth_output, file = expected_depth_output_path)

# save figures 
pdf(paste0(expected_depth_plots_prefix, "_boxplot.pdf"), width = 5)
print(expected_depth_all_snps$exp_dp_boxplot_plot)
dev.off()
