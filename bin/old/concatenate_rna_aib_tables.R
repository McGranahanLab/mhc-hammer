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

# get unique reads rna aib
single_mapping_aib_tables <- paste0(sample_name, ".", genes, ".", aligner, ".", snp_type, ".single_mapping.aib.csv")
all_single_mapping_aib_tables <- data.table()
for(single_mapping_aib_table in single_mapping_aib_tables){
  x <- fread(single_mapping_aib_table)
  all_single_mapping_aib_tables <- rbindlist(list(all_single_mapping_aib_tables, x))  
}
all_single_mapping_aib_tables[,patient := NULL]
all_single_mapping_aib_tables[,gene := gsub("(hla_.*?)_.*", "\\1", allele1)]
all_single_mapping_aib_tables[,gene := toupper(gene)]
all_single_mapping_aib_tables[,gene := gsub("_", "-", gene)]

setnames(all_single_mapping_aib_tables,
         old = c("paired_t_test", "paired_wilcoxon_test", "allele1_n_snps_with_coverage", "allele2_n_snps_with_coverage"),
         new = c("rna_aib_paired_t_test", "rna_aib_paired_wilcoxon_test", "rna_aib_allele1_n_snps_with_coverage", "rna_aib_allele2_n_snps_with_coverage"))

setcolorder(all_single_mapping_aib_tables, 
            c("sample_name", "gene", "allele1", "allele2", "rna_aib_paired_t_test",
              "rna_aib_paired_wilcoxon_test", "rna_aib_allele1_n_snps_with_coverage", "rna_aib_allele2_n_snps_with_coverage"))

fwrite(all_single_mapping_aib_tables, 
       file = paste0(sample_name, "_", snp_type, "_", aligner, "_uniquely_mapping_reads_rna_aib.csv"))

# get all reads rna aib
aib_tables <- paste0(sample_name, ".", genes, ".", aligner, ".", snp_type, ".aib.csv")
all_reads_aib_tables <- data.table()
for(aib_table in aib_tables){
  x <- fread(aib_table)
  all_reads_aib_tables <- rbindlist(list(all_reads_aib_tables, x))  
}
all_reads_aib_tables[,patient := NULL]
all_reads_aib_tables[,gene := gsub("(hla_.*?)_.*", "\\1", allele1)]
all_reads_aib_tables[,gene := toupper(gene)]
all_reads_aib_tables[,gene := gsub("_", "-", gene)]

setnames(all_reads_aib_tables,
         old = c("paired_t_test", "paired_wilcoxon_test", "allele1_n_snps_with_coverage", "allele2_n_snps_with_coverage"),
         new = c("rna_aib_paired_t_test", "rna_aib_paired_wilcoxon_test", "rna_aib_allele1_n_snps_with_coverage", "rna_aib_allele2_n_snps_with_coverage"))

setcolorder(all_reads_aib_tables, 
            c("sample_name", "gene", "allele1", "allele2", "rna_aib_paired_t_test",
              "rna_aib_paired_wilcoxon_test", "rna_aib_allele1_n_snps_with_coverage", "rna_aib_allele2_n_snps_with_coverage"))

fwrite(all_reads_aib_tables, 
       file = paste0(sample_name, "_", snp_type, "_", aligner, "_all_reads_rna_aib.csv"))

