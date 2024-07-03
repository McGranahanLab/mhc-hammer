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

# get the tables
cohort_splicing_tables <- data.table()
for(line_idx in 1:nrow(output_tables)){
  x <- fread(output_tables[line_idx]$csv_path)
  cohort_splicing_tables <- rbindlist(list(cohort_splicing_tables, x), use.names=TRUE, fill = TRUE)
}

# fix gene name
cohort_splicing_tables[,gene := toupper(gene)]
cohort_splicing_tables[,gene := gsub("_", "-", gene)]

# fix allele name
cohort_splicing_tables[,allele := toupper(allele)]
cohort_splicing_tables[,allele := gsub("HLA_", "HLA-", allele)]
cohort_splicing_tables[,allele := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele)]
cohort_splicing_tables[,allele := gsub("_", ":", allele)]

cohort_splicing_tables[, short_allele := gsub("(HLA-.*?:.*?):.*", "\\1", allele)]

# add consequence
cohort_splicing_tables[framshift == FALSE & premature_stop == TRUE, sj_consequence := "Inframe - PTC" ]
cohort_splicing_tables[framshift == TRUE & premature_stop == TRUE, sj_consequence := "Frameshift - PTC" ]
cohort_splicing_tables[framshift == TRUE & premature_stop == FALSE, sj_consequence := "Frameshift - no PTC" ]
cohort_splicing_tables[framshift == FALSE & premature_stop == FALSE, sj_consequence := "Inframe - no PTC" ]

# only keep some columns
cohort_splicing_tables <- cohort_splicing_tables[,c("tumour_sample_name", "normal_sample_name",
                                                    "gene", "allele", "short_allele",
                                                    "start", "end", "sj_consequence", "exon_intron_name", "sj_type",
                                                    "tumour_n_unique_reads", "normal_n_unique_reads",
                                                    "tumour_canonical_sj_read_count", "normal_canonical_sj_read_count",
                                                    "tumour_novel_transcript_proportion", "normal_novel_transcript_proportion",
                                                    "novel_transcript_proportion_change", "novel_transcript_proportion_changev2",
                                                    "fisher_pvalue", "fisher_odds_ratio", "fisher_ci_lower", "fisher_ci_upper")]

# add in the patient 
inventory <- fread(inventory_path)
rna_samples <- inventory[sequencing_type == "rnaseq", c("patient", "sample_name")]

cohort_splicing_tables <- merge(cohort_splicing_tables, rna_samples, 
                                by.x = "tumour_sample_name", by.y = "sample_name", all.x = TRUE)

if(nrow(cohort_splicing_tables[is.na(patient)]) > 0){
  stop("Why is the patient missing?")
}                     

# add in purity if it exists
tumour_wes_inventory <- inventory[sequencing_type == "wxs" & sample_type == "tumour",c("sample_name", "purity")]
if(nrow(tumour_wes_inventory) > 0){
  cohort_splicing_tables <- merge(cohort_splicing_tables, tumour_wes_inventory, 
                                  by.x = "tumour_sample_name", by.y = "sample_name", all.x = TRUE)
  cohort_splicing_tables[,purity_scaled_tumour_novel_transcript_proportion := tumour_novel_transcript_proportion/purity]
  cohort_splicing_tables[purity_scaled_tumour_novel_transcript_proportion > 1, purity_scaled_tumour_novel_transcript_proportion := 1]
  
  setcolorder(cohort_splicing_tables,
              c("patient", "tumour_sample_name", "normal_sample_name", 
                "gene", "allele", "short_allele", "sj_consequence",
                "start", "end", "tumour_n_unique_reads", "normal_n_unique_reads",
                "tumour_canonical_sj_read_count", "normal_canonical_sj_read_count",
                "tumour_novel_transcript_proportion", "normal_novel_transcript_proportion",
                "novel_transcript_proportion_change", "novel_transcript_proportion_changev2",
                "fisher_pvalue", "fisher_odds_ratio", "fisher_ci_lower", "fisher_ci_upper",
                "purity", "purity_scaled_tumour_novel_transcript_proportion"))
  
}else{
  setcolorder(cohort_splicing_tables,
              c("patient", "tumour_sample_name", "normal_sample_name", 
                "gene", "allele", "short_allele", "sj_consequence",
                "start", "end", "tumour_n_unique_reads", "normal_n_unique_reads",
                "tumour_canonical_sj_read_count", "normal_canonical_sj_read_count",
                "tumour_novel_transcript_proportion", "normal_novel_transcript_proportion",
                "novel_transcript_proportion_change", "novel_transcript_proportion_changev2",
                "fisher_pvalue", "fisher_odds_ratio", "fisher_ci_lower", "fisher_ci_upper"))
}

fwrite(cohort_splicing_tables, "tumour_normal_splicing_events.csv")
