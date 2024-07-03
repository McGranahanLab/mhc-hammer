library(data.table)
library(argparse)

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--alleles',
                    help='Genes with output',
                    required=TRUE, nargs = "+")
parser$add_argument('--snp_type',
                    help="Either 'exon_snps' or 'all_snps'",
                    required=TRUE)
parser$add_argument('--tumour_sample_name',
                    help="Tumour sample name'",
                    required=TRUE)
parser$add_argument('--normal_sample_name',
                    help="Normal sample name'",
                    required=TRUE)
parser$add_argument('--aligner',
                    help="novoalign for DNA",
                    required=TRUE)

args <- parser$parse_args()
alleles <- args$alleles
snp_type <- args$snp_type
tumour_sample_name <- args$tumour_sample_name
normal_sample_name <- args$normal_sample_name
aligner <- args$aligner

cat("alleles=", alleles, "\n")
cat("snp_type=", snp_type, "\n")
cat("tumour_sample_name=", tumour_sample_name, "\n")
cat("normal_sample_name=", normal_sample_name, "\n")
cat("aligner=", aligner, "\n")

repression_tables <- paste0(tumour_sample_name, ".", alleles, ".", aligner, ".", snp_type, ".tumour_normal_comparison.csv")
all_repression_tables <- data.table()
for(repression_table in repression_tables){
  x <- fread(repression_table)
  all_repression_tables <- rbindlist(list(all_repression_tables, x))  
}

all_repression_tables[, patient := NULL]
all_repression_tables[, normal_sample_name := normal_sample_name]

setnames(all_repression_tables, 
         c("paired_t_test", "paired_wilcoxon_test", "median_tumour_dp", "median_normal_dp", "normal_n_snps_with_coverage", "tumour_n_snps_with_coverage"),
         c("repression_paired_t_test", "repression_paired_wilcoxon_test", "repression_median_tumour_dp", "repression_median_normal_dp", "repression_normal_n_snps_with_coverage", "repression_tumour_n_snps_with_coverage"))

fwrite(all_repression_tables, file = paste0(tumour_sample_name, "_", snp_type, "_", aligner, "_rna_repression.csv"))

