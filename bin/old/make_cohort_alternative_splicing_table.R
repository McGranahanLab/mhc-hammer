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

output_tables[grepl("known", csv_path), type := "known_sjs"]
output_tables[grepl("novel", csv_path), type := "novel_sjs"]
output_tables[grepl("tumour_normal", csv_path), type := "tumour_normal_sjs"]

# check if tumour has matched normal
inventory <- fread(inventory_path)
novel_sjs_paths <- inventory[sequencing_type == "rnaseq",c("sample_name", "normal_sample_name")]
novel_sjs_paths[normal_sample_name == "", novel_sj_path := paste0(sample_name, "_novel_splice_junctions.csv")]
novel_sjs_paths[normal_sample_name != "", novel_sj_path := paste0(sample_name, "_tumour_normal_splice_junctions.csv")]

# get the novel tables
cohort_novel_splicing_tables <- data.table()
for(line_idx in 1:nrow(novel_sjs_paths)){
  
  cat(line_idx, "/", nrow(novel_sjs_paths), "\n")
  
  x <- fread(novel_sjs_paths[line_idx]$novel_sj_path)
  
  if(novel_sjs_paths[line_idx]$normal_sample_name != ""){
    x[,sample_name := gsub("_tumour_normal_splice_junctions.csv$", "", basename(novel_sjs_paths[line_idx]$novel_sj_path))]
    x[,matched_normal_sample_name := novel_sjs_paths[line_idx]$normal_sample_name]
    setnames(x, 
             old = c("tumour_n_unique_reads", "tumour_n_multtimap_reads", "tumour_max_overhang", "tumour_intron_n_reads", "tumour_canonical_sj_read_count", "tumour_total_read_count", "tumour_novel_transcript_proportion"),
             new = c("n_unique_reads", "n_multtimap_reads", "max_overhang", "intron_n_reads", "canonical_sj_read_count", "total_read_count", "novel_transcript_proportion"))
  }else{
    x[,sample_name := gsub("_novel_splice_junctions.csv$", "", basename(novel_sjs_paths[line_idx]$novel_sj_path))]
  }
  
  cohort_novel_splicing_tables <- rbindlist(list(cohort_novel_splicing_tables, x), use.names=TRUE, fill = TRUE)
  
}

# get the known tables
known_tables <- output_tables[type == "known_sjs"]
cohort_known_splicing_tables <- data.table()
for(line_idx in 1:nrow(known_tables)){
  
  cat(line_idx, "/", nrow(known_tables), "\n")
  
  x <- fread(known_tables[line_idx]$csv_path)
  
  x[,sample_name := gsub("_known_splice_junctions.csv$", "", basename(known_tables[line_idx]$csv_path))]
  
  cohort_known_splicing_tables <- rbindlist(list(cohort_known_splicing_tables, x), use.names=TRUE)
  
}

# update the gene name
cohort_novel_splicing_tables[,gene := toupper(gene)]
cohort_novel_splicing_tables[,gene := gsub("_", "-", gene)]
cohort_known_splicing_tables[,gene := toupper(gene)]
cohort_known_splicing_tables[,gene := gsub("_", "-", gene)]

# update the allele name
cohort_novel_splicing_tables[,allele := toupper(allele)]
cohort_novel_splicing_tables[,allele := gsub("HLA_", "HLA-", allele)]
cohort_novel_splicing_tables[,allele := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele)]
cohort_novel_splicing_tables[,allele := gsub("_", ":", allele)]

cohort_known_splicing_tables[,allele := toupper(allele)]
cohort_known_splicing_tables[,allele := gsub("HLA_", "HLA-", allele)]
cohort_known_splicing_tables[,allele := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele)]
cohort_known_splicing_tables[,allele := gsub("_", ":", allele)]

# add in the short allele
cohort_novel_splicing_tables[, short_allele := gsub("(HLA-.*?:.*?):.*", "\\1", allele)]
cohort_known_splicing_tables[, short_allele := gsub("(HLA-.*?:.*?):.*", "\\1", allele)]

# add consequence
cohort_novel_splicing_tables[framshift == FALSE & premature_stop == TRUE, sj_consequence := "Inframe - PTC" ]
cohort_novel_splicing_tables[framshift == TRUE & premature_stop == TRUE, sj_consequence := "Frameshift - PTC" ]
cohort_novel_splicing_tables[framshift == TRUE & premature_stop == FALSE, sj_consequence := "Frameshift - no PTC" ]
cohort_novel_splicing_tables[framshift == FALSE & premature_stop == FALSE, sj_consequence := "Inframe" ]

# only keep some columns
cohort_novel_splicing_tables <- cohort_novel_splicing_tables[,c("sample_name", "matched_normal_sample_name",
                                                                "gene", "allele", "short_allele",
                                                                "start", "end", 
                                                                "n_unique_reads", "normal_n_unique_reads",
                                                                "sj_in_feature",
                                                                "canonical_sj_read_count", 
                                                                "normal_canonical_sj_read_count",
                                                                "total_read_count", "normal_total_read_count",
                                                                "novel_transcript_proportion", "normal_novel_transcript_proportion" ,
                                                                "novel_transcript_proportion_change", "novel_transcript_proportion_changev2",
                                                                "fisher_pvalue", "fisher_odds_ratio",
                                                                "fisher_ci_lower", "fisher_ci_upper", 
                                                                "sj_type", "exon_intron_name", 
                                                                "premature_stop", "framshift",
                                                                "sj_consequence")]

cohort_known_splicing_tables <- cohort_known_splicing_tables[,c("sample_name", "allele", "short_allele", "gene",
                                                                "start", "end", "n_unique_reads", "known_feature_name")]

# add in the patient and sample_type
rna_samples <- inventory[sequencing_type == "rnaseq",c("patient", "sample_name", "sample_type")]
cohort_novel_splicing_tables <- merge(cohort_novel_splicing_tables, rna_samples, 
                                      by = "sample_name", all.x = TRUE)
cohort_known_splicing_tables <- merge(cohort_known_splicing_tables, rna_samples, 
                                      by = "sample_name", all.x = TRUE)

if(nrow(cohort_novel_splicing_tables[is.na(patient)]) > 0){
  stop("Why is the patient missing?")
}

if(nrow(cohort_known_splicing_tables[is.na(patient)]) > 0){
  stop("Why is the patient missing?")
}

# add in if sample passes in rna and wes
# allele_table <- fread(allele_overview_path)
# allele_table <- unique(allele_table[,c("sample_name", "allele", "rna_fail", "wes_fail")])
# cohort_novel_splicing_tables <- merge(cohort_novel_splicing_tables, allele_table, by = c("sample_name", "allele"), all.x = TRUE)

# add in the purity and the purity scaled novel_transcript_proportion
tumour_wes_inventory <- inventory[sequencing_type == "wxs" & sample_type == "tumour",c("sample_name", "purity")]
if(nrow(tumour_wes_inventory) > 0){
  cohort_novel_splicing_tables <- merge(cohort_novel_splicing_tables, tumour_wes_inventory, 
                                        by = "sample_name", all.x = TRUE)
  cohort_known_splicing_tables <- merge(cohort_known_splicing_tables, tumour_wes_inventory, 
                                        by = "sample_name", all.x = TRUE)
  cohort_novel_splicing_tables[,purity_scaled_novel_transcript_proportion := novel_transcript_proportion/purity]
  cohort_novel_splicing_tables[purity_scaled_novel_transcript_proportion > 1, purity_scaled_novel_transcript_proportion := 1]
  
}

setcolorder(cohort_novel_splicing_tables,
            c("patient", "sample_name", "matched_normal_sample_name", "gene", "allele", "short_allele", "start", "end", "n_unique_reads", "normal_n_unique_reads", "canonical_sj_read_count", "normal_canonical_sj_read_count", "total_read_count", "normal_total_read_count", "novel_transcript_proportion", "normal_novel_transcript_proportion", "novel_transcript_proportion_change", "novel_transcript_proportion_changev2", "fisher_pvalue", "fisher_odds_ratio", "fisher_ci_lower", "fisher_ci_upper", "sj_type", "exon_intron_name", "premature_stop", "framshift", "sj_consequence", "sample_type", "purity", "purity_scaled_novel_transcript_proportion"))

setcolorder(cohort_known_splicing_tables,
            c("patient", "sample_name", "gene", "allele", "short_allele", "start", "end", "n_unique_reads", "known_feature_name", "sample_type", "purity"))

fwrite(cohort_novel_splicing_tables, "novel_splicing_events.csv")
fwrite(cohort_known_splicing_tables, "known_splicing_events.csv")
