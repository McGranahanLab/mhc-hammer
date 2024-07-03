library(data.table)
library(argparse)

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--inventory_path',
                    help='Path to inventory.',
                    required=TRUE)
parser$add_argument('--csv_tables_path',
                    help="Path to all the patient's *rpkm.csv",
                    required=FALSE)

args <- parser$parse_args()
inventory_path <- args$inventory_path
csv_tables_path <- args$csv_tables_path

cat("inventory_path=", inventory_path, "\n")
cat("csv_tables_path=", csv_tables_path, "\n")

#### get output tables ####

output_tables <- fread(csv_tables_path, header = FALSE, col.names = "csv_path")
output_tables <- output_tables[grepl("known", csv_path)]
              
# get the known tables
cohort_known_splicing_tables <- data.table()
for(line_idx in 1:nrow(output_tables)){
  x <- fread(output_tables[line_idx]$csv_path)
  cohort_known_splicing_tables <- rbindlist(list(cohort_known_splicing_tables, x), use.names=TRUE)
}

# update the gene name
cohort_known_splicing_tables[,gene := toupper(gene)]
cohort_known_splicing_tables[,gene := gsub("_", "-", gene)]

# update the allele name
cohort_known_splicing_tables[,allele := toupper(allele)]
cohort_known_splicing_tables[,allele := gsub("HLA_", "HLA-", allele)]
cohort_known_splicing_tables[,allele := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele)]
cohort_known_splicing_tables[,allele := gsub("_", ":", allele)]

# add in the short allele
cohort_known_splicing_tables[, short_allele := gsub("(HLA-.*?:.*?):.*", "\\1", allele)]

# remove the columns we dont need
cohort_known_splicing_tables <- cohort_known_splicing_tables[,c("sample_name", "allele", "short_allele", "gene",
                                                                "start", "end", "n_unique_reads", "known_feature_name")]

# add in the patient and sample_type from the inventory
inventory <- fread(inventory_path)
rna_samples <- inventory[sequencing_type == "rnaseq",c("patient", "sample_name", "sample_type")]
cohort_known_splicing_tables <- merge(cohort_known_splicing_tables, rna_samples, 
                                      by = "sample_name", all.x = TRUE)

if(nrow(cohort_known_splicing_tables[is.na(patient)]) > 0){
  stop("Why is the patient missing?")
}

setcolorder(cohort_known_splicing_tables,
            c("patient", "sample_name", "gene", "allele", "short_allele", "start", 
              "end", "n_unique_reads", "known_feature_name", "sample_type"))

fwrite(cohort_known_splicing_tables, "known_splicing_events.csv")
