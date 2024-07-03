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

output_tables[grepl("wxs.library_size_without_unmapped.txt$", csv_path), table_type := "wes_without_unmapped"]
output_tables[grepl("wxs.library_size_with_unmapped.txt$", csv_path), table_type := "wes_with_unmapped"]
output_tables[grepl("rnaseq.library_size_without_unmapped.txt$", csv_path), table_type := "rna_without_unmapped"]
output_tables[grepl("rnaseq.library_size_with_unmapped.txt$", csv_path), table_type := "rna_with_unmapped"]

wes_without_unmapped_tables <- output_tables[table_type == "wes_without_unmapped"]
wes_with_unmapped_tables <- output_tables[table_type == "wes_with_unmapped"]
rna_without_unmapped_tables <- output_tables[table_type == "rna_without_unmapped"]
rna_with_unmapped_tables <- output_tables[table_type == "rna_with_unmapped"]

# wes without unmapped
cohort_wes_without_unmapped <- data.table()
for(line_idx in 1:nrow(wes_without_unmapped_tables)){
  
  x <- as.numeric(readLines(wes_without_unmapped_tables[line_idx]$csv_path))
  
  cohort_wes_without_unmapped <- rbindlist(list(cohort_wes_without_unmapped, 
                                                data.table(csv_path = wes_without_unmapped_tables[line_idx]$csv_path,
                                                           read_count_without_unmapped = x)))
  
}

# wes with unmapped
cohort_wes_with_unmapped <- data.table()
for(line_idx in 1:nrow(wes_with_unmapped_tables)){
  
  x <- as.numeric(readLines(wes_with_unmapped_tables[line_idx]$csv_path))
  
  cohort_wes_with_unmapped <- rbindlist(list(cohort_wes_with_unmapped, 
                                             data.table(csv_path = wes_with_unmapped_tables[line_idx]$csv_path,
                                                        read_count_with_unmapped = x)))
  
}

# rna without unmapped
cohort_rna_without_unmapped <- data.table()
for(line_idx in 1:nrow(rna_without_unmapped_tables)){
  
  x <- as.numeric(readLines(rna_without_unmapped_tables[line_idx]$csv_path))
  
  cohort_rna_without_unmapped <- rbindlist(list(cohort_rna_without_unmapped, 
                                                data.table(csv_path = rna_without_unmapped_tables[line_idx]$csv_path,
                                                           read_count_without_unmapped = x)))
  
}

# rna with unmapped
cohort_rna_with_unmapped <- data.table()
for(line_idx in 1:nrow(rna_with_unmapped_tables)){
  
  x <- as.numeric(readLines(rna_with_unmapped_tables[line_idx]$csv_path))
  
  cohort_rna_with_unmapped <- rbindlist(list(cohort_rna_with_unmapped, 
                                             data.table(csv_path = rna_with_unmapped_tables[line_idx]$csv_path,
                                                        read_count_with_unmapped = x)))
  
}

cohort_wes_with_unmapped[, sample_name := gsub("_wxs.library_size_with_unmapped.txt", "", csv_path)]
cohort_wes_with_unmapped[,csv_path := NULL]

cohort_wes_without_unmapped[, sample_name := gsub("_wxs.library_size_without_unmapped.txt", "", csv_path)]
cohort_wes_without_unmapped[,csv_path := NULL]

cohort_rna_with_unmapped[, sample_name := gsub("_rnaseq.library_size_with_unmapped.txt", "", csv_path)]
cohort_rna_with_unmapped[,csv_path := NULL]

cohort_rna_without_unmapped[, sample_name := gsub("_rnaseq.library_size_without_unmapped.txt", "", csv_path)]
cohort_rna_without_unmapped[,csv_path := NULL]

wes_library_size_dt <- merge(cohort_wes_with_unmapped, cohort_wes_without_unmapped, by = "sample_name", all = TRUE)
if(nrow(wes_library_size_dt[is.na(read_count_without_unmapped) | is.na(read_count_with_unmapped)]) > 0){
  stop("why are wes samples missing a library size")
}
wes_library_size_dt[, fraction_align := read_count_without_unmapped/read_count_with_unmapped]

rna_library_size_dt <- merge(cohort_rna_with_unmapped, cohort_rna_without_unmapped, by = "sample_name", all = TRUE)
if(nrow(rna_library_size_dt[is.na(read_count_without_unmapped) | is.na(read_count_with_unmapped)]) > 0){
  stop("why are wes samples missing a library size")
}
rna_library_size_dt[, fraction_align := read_count_without_unmapped/read_count_with_unmapped]

wes_library_size_dt[,sequencing_type := "wes"]
rna_library_size_dt[,sequencing_type := "rnaseq"]

library_size_dt <- rbindlist(list(wes_library_size_dt, rna_library_size_dt))

fwrite(library_size_dt, "cohort_library_size.csv")
