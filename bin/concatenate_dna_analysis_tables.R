library(data.table)
library(argparse)

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--genes',
                    help='Genes with output',
                    required=TRUE, nargs = "+")
parser$add_argument('--snp_type',
                    help="Either 'exon_snps' or 'all_snps'",
                    required=TRUE)
parser$add_argument('--sample_name',
                    help="Sample name'",
                    required=TRUE)
parser$add_argument('--aligner',
                    help="novoalign for DNA",
                    required=TRUE)

args <- parser$parse_args()
genes <- args$genes
snp_type <- args$snp_type
sample_name <- args$sample_name
aligner <- args$aligner

cat("genes=", genes, "\n")
cat("snp_type=", snp_type, "\n")
cat("sample_name=", sample_name, "\n")
cat("aligner=", aligner, "\n")

cn_tables <- paste0(sample_name, ".", genes, ".", snp_type, ".cn.csv")
all_cn_tables <- data.table()
for(cn_table in cn_tables){
  x <- fread(cn_table)
  all_cn_tables <- rbindlist(list(all_cn_tables, x))  
}

aib_tables <- paste0(sample_name, ".", genes, ".", snp_type, ".logr_aib.csv")
all_aib_tables <- data.table()
for(aib_table in aib_tables){
  x <- fread(aib_table)
  all_aib_tables <- rbindlist(list(all_aib_tables, x))  
}

expected_depth_tables <- paste0(sample_name, ".", genes, ".", snp_type, ".expected_depth.csv")
all_expected_depth_tables <- data.table()
for(expected_depth_table in expected_depth_tables){
  x <- fread(expected_depth_table)
  all_expected_depth_tables <- rbindlist(list(all_expected_depth_tables, x))  
}

# merge the tables together
all_expected_depth_tables[, gene := gsub("(hla_.*?)_.*", "\\1", allele1)]
all_expected_depth_tables[,c("allele1", "allele2") := NULL]

all_aib_tables[, gene := gsub("(hla_.*?)_.*", "\\1", allele1)]
all_aib_tables[,c("allele1", "allele2") := NULL]

all_cn_tables[, gene := gsub("(hla_.*?)_.*", "\\1", allele1)]

setnames(all_cn_tables, c("n_snps", "n_bins"), c("cn_n_snps", "cn_n_bins"))
setnames(all_aib_tables, 
         c("n_snps", "paired_t_test", "paired_wilcoxon_test"), 
         c("logr_aib_n_snps", "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test"))

dna_analysis_dt <- merge(all_expected_depth_tables, all_aib_tables, by = "gene", all = TRUE)
dna_analysis_dt <- merge(dna_analysis_dt, all_cn_tables, by = "gene", all = TRUE)
dna_analysis_dt[,sample_name := sample_name]
setcolorder(dna_analysis_dt, 
            c("sample_name", "gene", "allele1", "allele2",
              "allele1_exp_dp", "allele2_exp_dp",
              "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test", "logr_aib_n_snps",
              "cn1", "cn1_lower", "cn1_upper",
              "cn2", "cn2_lower", "cn2_upper",
              "cn1_binned", "cn1_binned_lower", "cn1_binned_upper",
              "cn2_binned", "cn2_binned_lower", "cn2_binned_upper",
              "cn_n_bins", "cn_n_snps"))

dna_analysis_dt[,gene := toupper(gene)]
dna_analysis_dt[,gene := gsub("_", "-", gene)]

fwrite(dna_analysis_dt, file = paste0(sample_name, "_", snp_type, "_", aligner, "_dna_analysis.csv"))

