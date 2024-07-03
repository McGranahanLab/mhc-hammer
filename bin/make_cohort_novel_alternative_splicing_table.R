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
output_tables <- output_tables[grepl("novel", csv_path)]

# get the novel tables
cohort_novel_splicing_tables <- data.table()
for(line_idx in 1:nrow(output_tables)){
  x <- fread(output_tables[line_idx]$csv_path)
  cohort_novel_splicing_tables <- rbindlist(list(cohort_novel_splicing_tables, x), use.names=TRUE, fill = TRUE)
}

# fix the gene name
cohort_novel_splicing_tables[,gene := toupper(gene)]
cohort_novel_splicing_tables[,gene := gsub("_", "-", gene)]

# fix the allele name
cohort_novel_splicing_tables[,allele := toupper(allele)]
cohort_novel_splicing_tables[,allele := gsub("HLA_", "HLA-", allele)]
cohort_novel_splicing_tables[,allele := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele)]
cohort_novel_splicing_tables[,allele := gsub("_", ":", allele)]

cohort_novel_splicing_tables[, short_allele := gsub("(HLA-.*?:.*?):.*", "\\1", allele)]

# add consequence
cohort_novel_splicing_tables[framshift == FALSE & premature_stop == TRUE, sj_consequence := "Inframe - PTC" ]
cohort_novel_splicing_tables[framshift == TRUE & premature_stop == TRUE, sj_consequence := "Frameshift - PTC" ]
cohort_novel_splicing_tables[framshift == TRUE & premature_stop == FALSE, sj_consequence := "Frameshift - no PTC" ]
cohort_novel_splicing_tables[framshift == FALSE & premature_stop == FALSE, sj_consequence := "Inframe - no PTC" ]

# only keep some columns
cohort_novel_splicing_tables <- cohort_novel_splicing_tables[,c("sample_name", 
                                                                "gene", "allele", "short_allele",
                                                                "start", "end", 
                                                                "n_unique_reads", 
                                                                "sj_in_feature",
                                                                "canonical_sj_read_count", 
                                                                "total_read_count", 
                                                                "novel_transcript_proportion", 
                                                                "sj_type", "exon_intron_name", 
                                                                "premature_stop", "framshift",
                                                                "sj_consequence")]

# add in the patient and sample type
inventory <- fread(inventory_path)
rna_samples <- inventory[sequencing_type == "rnaseq",c("patient", "sample_name", "sample_type")]

cohort_novel_splicing_tables <- merge(cohort_novel_splicing_tables, rna_samples, 
                                      by = "sample_name", all.x = TRUE)

if(nrow(cohort_novel_splicing_tables[is.na(patient)]) > 0){
  stop("Why is the patient missing?")
}                     

# add in the purity 
tumour_wes_inventory <- inventory[sequencing_type == "wxs" & sample_type == "tumour",c("sample_name", "purity")]
if(nrow(tumour_wes_inventory) > 0){
  cohort_novel_splicing_tables <- merge(cohort_novel_splicing_tables, tumour_wes_inventory, 
                                        by = "sample_name", all.x = TRUE)
  cohort_novel_splicing_tables[,purity_scaled_novel_transcript_proportion := novel_transcript_proportion/purity]
  cohort_novel_splicing_tables[purity_scaled_novel_transcript_proportion > 1, purity_scaled_novel_transcript_proportion := 1]
  
  setcolorder(cohort_novel_splicing_tables,
              c("patient", "sample_name", "gene", "allele", "short_allele", 
                "start", "end", "n_unique_reads", "canonical_sj_read_count", 
                "total_read_count", "novel_transcript_proportion",
                "sj_type", "exon_intron_name", "premature_stop", "framshift", 
                "sj_consequence", "sample_type", "purity", "purity_scaled_novel_transcript_proportion"))
  
}else{
  setcolorder(cohort_novel_splicing_tables,
              c("patient", "sample_name", "gene", "allele", "short_allele", 
                "start", "end", "n_unique_reads", "canonical_sj_read_count", 
                "total_read_count", "novel_transcript_proportion",
                "sj_type", "exon_intron_name", "premature_stop", "framshift", 
                "sj_consequence", "sample_type"))
}


fwrite(cohort_novel_splicing_tables, "novel_splicing_events.csv")

