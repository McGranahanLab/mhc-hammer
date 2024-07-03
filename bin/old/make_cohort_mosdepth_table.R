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

# mosdepth tables
output_tables[grepl("rnaseq.star.mosdepth.bed$", csv_path), table_type := "rnaseq_star"]
output_tables[grepl("rnaseq.novoalign.mosdepth.bed$", csv_path), table_type := "rnaseq_novoalign"]
output_tables[grepl("wxs.novoalign.mosdepth.bed$", csv_path), table_type := "wxs_novoalign"]

rnaseq_star_mosdepth_tables <- output_tables[table_type == "rnaseq_star"]
rnaseq_novoalign_mosdepth_tables <- output_tables[table_type == "rnaseq_novoalign"]
wxs_novoalign_mosdepth_tables <- output_tables[table_type == "wxs_novoalign"]

# rna star
cohort_rnaseq_star_mosdepth_tables <- data.table()
for(line_idx in 1:nrow(rnaseq_star_mosdepth_tables)){
  
  x <- fread(rnaseq_star_mosdepth_tables[line_idx]$csv_path)
  setnames(x, c("allele", "start", "stop", "feature_name", "depth"))
  
  x[,sample_name := gsub(".rnaseq.star.mosdepth.bed", "", basename(rnaseq_star_mosdepth_tables[line_idx]$csv_path))]
  
  cohort_rnaseq_star_mosdepth_tables <- rbindlist(list(cohort_rnaseq_star_mosdepth_tables, x))
  
}

# rna novoalign
cohort_rnaseq_novoalign_mosdepth_tables <- data.table()
for(line_idx in 1:nrow(rnaseq_novoalign_mosdepth_tables)){
  
  x <- fread(rnaseq_novoalign_mosdepth_tables[line_idx]$csv_path)
  setnames(x, c("allele", "start", "stop", "feature_name", "depth"))
  
  x[,sample_name := gsub(".rnaseq.novoalign.mosdepth.bed", "", basename(rnaseq_novoalign_mosdepth_tables[line_idx]$csv_path))]
  
  cohort_rnaseq_novoalign_mosdepth_tables <- rbindlist(list(cohort_rnaseq_novoalign_mosdepth_tables, x))
  
}

# wxs novoalign
cohort_wxs_novoalign_mosdepth_tables <- data.table()
for(line_idx in 1:nrow(wxs_novoalign_mosdepth_tables)){
  
  x <- fread(wxs_novoalign_mosdepth_tables[line_idx]$csv_path)
  setnames(x, c("allele", "start", "stop", "feature_name", "depth"))
  
  x[,sample_name := gsub(".wxs.novoalign.mosdepth.bed", "", basename(wxs_novoalign_mosdepth_tables[line_idx]$csv_path))]
  
  cohort_wxs_novoalign_mosdepth_tables <- rbindlist(list(cohort_wxs_novoalign_mosdepth_tables, x))
  
}

setnames(cohort_rnaseq_star_mosdepth_tables, "stop", "end")
setnames(cohort_rnaseq_novoalign_mosdepth_tables, "stop", "end")
setnames(cohort_wxs_novoalign_mosdepth_tables, "stop", "end")

setcolorder(cohort_rnaseq_star_mosdepth_tables, c("sample_name", "allele", "start", "end", "feature_name", "depth"))
setcolorder(cohort_rnaseq_novoalign_mosdepth_tables, c("sample_name", "allele", "start", "end", "feature_name", "depth"))
setcolorder(cohort_wxs_novoalign_mosdepth_tables, c("sample_name", "allele", "start", "end", "feature_name", "depth"))

fwrite(cohort_rnaseq_star_mosdepth_tables, "mosdepth_rnaseq_star.csv")
fwrite(cohort_rnaseq_novoalign_mosdepth_tables, "mosdepth_novoalign_rnaseq_star.csv")
fwrite(cohort_wxs_novoalign_mosdepth_tables, "mosdepth_novoalign_wes_star.csv")
