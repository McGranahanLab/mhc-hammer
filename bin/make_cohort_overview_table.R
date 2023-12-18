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
parser$add_argument('--hlahd_germline_samples_path',
                    help="Path to the germline samples used by HLAHD",
                    required=FALSE)
parser$add_argument('--max_cn_range',
                    help="Maximum range in cn",
                    required=TRUE)
parser$add_argument('--min_n_snps',
                    help="Minimum number of snps",
                    required=TRUE)
parser$add_argument('--min_expected_depth',
                    help="Minimum number of snps",
                    required=TRUE)

parser$add_argument('--min_frac_mapping_uniquely',
                    help="",
                    required=TRUE)
parser$add_argument('--max_frac_mapping_multi_gene',
                    help="",
                    required=TRUE)

parser$add_argument('--dna_snp_min_depth',
                    help="",
                    required=TRUE)

args <- parser$parse_args()
inventory_path <- args$inventory_path
csv_tables_path <- args$csv_tables_path
hlahd_germline_samples_path <- args$hlahd_germline_samples_path
outfile <- args$outfile
max_cn_range <- args$max_cn_range
min_n_snps <- args$min_n_snps
min_expected_depth <- args$min_expected_depth
min_frac_mapping_uniquely <- args$min_frac_mapping_uniquely
max_frac_mapping_multi_gene <- args$max_frac_mapping_multi_gene
dna_snp_min_depth <- args$dna_snp_min_depth

cat("inventory_path=", inventory_path, "\n")
cat("csv_tables_path=", csv_tables_path, "\n")
cat("hlahd_germline_samples_path=", hlahd_germline_samples_path, "\n")
cat("dna_snp_min_depth=", dna_snp_min_depth, "\n")

#### get output tables ####

output_tables <- fread(csv_tables_path, header = FALSE, col.names = "csv_path")

# dna tables
output_tables[grepl("all_snps_novoalign_dna_analysis.csv$", csv_path), table_type := "dna_analysis"]
output_tables[grepl("genome_allele_table.csv$", csv_path), table_type := "genome_allele_table"]
output_tables[grepl("wxs_novoalign.hla_bam_read_count.csv$", csv_path), table_type := "dna_bam_read_count"]

# rna tables
output_tables[grepl("novoalign_rpkm.csv$", csv_path), table_type := "rpkm"]
output_tables[grepl("all_snps_novoalign_rna_repression.csv$", csv_path), table_type := "repression"]
output_tables[grepl("all_snps_novoalign_all_reads_rna_aib.csv$", csv_path), table_type := "rna_aib"]
output_tables[grepl("transcriptome_allele_table.csv$", csv_path), table_type := "transcriptome_allele_table"]
output_tables[grepl("_rnaseq_novoalign.hla_bam_read_count.csv$", csv_path), table_type := "rna_bam_read_count"]

# library_size
output_tables[csv_path == "cohort_library_size.csv", table_type := "library_size"]

#### Make the initial table with the samples with rna and dna ####
inventory <- fread(inventory_path)
setnames(inventory, "normal_sample_id", "normal_sample_name")

wes_overview_dt <- unique(inventory[sequencing_type == "wxs" & sample_type == "tumour",
                                    c("patient", "sample_name", "purity", "ploidy", "normal_sample_name", "sample_type")])
setnames(wes_overview_dt, "normal_sample_name", "wes_germline_id")
wes_overview_dt[, sample_has_wes := TRUE]

rna_overview_dt <- unique(inventory[sequencing_type == "rnaseq",
                                    c("patient", "sample_name", "normal_sample_name", "sample_type")])
setnames(rna_overview_dt, "normal_sample_name", "rnaseq_normal_sample_name")
rna_overview_dt[, sample_has_rna := TRUE]

overview_table <- merge(wes_overview_dt, rna_overview_dt, 
                        by = c("patient", "sample_name"),
                        all = TRUE)
overview_table[is.na(sample_type.x), sample_type.x := sample_type.y]
overview_table[,sample_type.y := NULL]
setnames(overview_table, "sample_type.x", "sample_type")

overview_table[is.na(sample_has_wes), sample_has_wes := FALSE]
overview_table[is.na(sample_has_rna), sample_has_rna := FALSE]

# remove tumour samples that only have rna
overview_table <- overview_table[!(sample_type == "tumour" & sample_has_wes == FALSE)]
overview_table[,rnaseq_normal_sample_name := as.character(rnaseq_normal_sample_name)]
overview_table[is.na(rnaseq_normal_sample_name), rnaseq_normal_sample_name := ""]

##### Add in the germline sample used for hlahd #####
hlahd_gl_samples <- fread(hlahd_germline_samples_path, header = FALSE)
hlahd_gl_sample_mapping <- unique(inventory[sample_name %in% hlahd_gl_samples$V1, c("patient", "sample_name")])

# check patient is there once
if(nrow(hlahd_gl_sample_mapping[,.N,by = "patient"][N>1]) > 0){
  stop("Patients should only have one hlahd output")
}

setnames(hlahd_gl_sample_mapping, "sample_name", "hlahd_germline_sample_name")
overview_table <- merge(overview_table, hlahd_gl_sample_mapping, by = "patient", all.x = TRUE)
if(nrow(overview_table[is.na(hlahd_germline_sample_name)]) > 0){
  stop("Patients are missing HLAHD?")
}

##### add in the library size #####
if(length(output_tables[table_type == "library_size"]$csv_path) != 1){
  stop("why no library size?")
}

library_size <- fread(output_tables[table_type == "library_size"]$csv_path)

# merge in the wes tumour
wes_library_size <- library_size[sequencing_type == "wes",
                                 c("sample_name", "read_count_with_unmapped",
                                   "read_count_without_unmapped", "fraction_align")]
setnames(wes_library_size, 
         c("read_count_with_unmapped", "read_count_without_unmapped", "fraction_align"),
         c("wes_tumour_read_count_with_unmapped", "wes_tumour_read_count_without_unmapped", "wes_tumour_fraction_align"))

overview_table <- merge(overview_table, wes_library_size, 
                        by = c("sample_name"), all.x = TRUE)

# merge the wes germline
setnames(wes_library_size, 
         c("wes_tumour_read_count_with_unmapped", "wes_tumour_read_count_without_unmapped", "wes_tumour_fraction_align"),
         c("wes_germline_read_count_with_unmapped", "wes_germline_read_count_without_unmapped", "wes_germline_fraction_align"))

overview_table <- merge(overview_table, wes_library_size, 
                        by.x = "wes_germline_id", 
                        by.y = "sample_name", all.x = TRUE)

# merge the wes hlahd sample name
setnames(wes_library_size, 
         c("wes_germline_read_count_with_unmapped", "wes_germline_read_count_without_unmapped", "wes_germline_fraction_align"),
         c("wes_hlahd_germline_read_count_with_unmapped", "wes_hlahd_germline_read_count_without_unmapped", "wes_hlahd_germline_fraction_align"))

overview_table <- merge(overview_table, wes_library_size, 
                        by.x = "hlahd_germline_sample_name", 
                        by.y = "sample_name", all.x = TRUE)

# merge in the rnaseq sample
rna_library_size <- library_size[sequencing_type == "rnaseq",
                                 c("sample_name", "read_count_with_unmapped",
                                   "read_count_without_unmapped", "fraction_align")]
setnames(rna_library_size, 
         c("read_count_with_unmapped", "read_count_without_unmapped", "fraction_align"),
         c("rna_read_count_with_unmapped", "rna_read_count_without_unmapped", "rna_readcount_fraction_align"))

overview_table <- merge(overview_table, rna_library_size, 
                        by = c("sample_name"), all.x = TRUE)

# merge in matched normal 
setnames(rna_library_size, 
         c("rna_read_count_with_unmapped", "rna_read_count_without_unmapped", "rna_readcount_fraction_align"),
         c("rna_normal_read_count_with_unmapped", "rna_normal_read_count_without_unmapped", "rna_normal_readcount_fraction_align"))

overview_table <- merge(overview_table, rna_library_size, 
                        by.x = c("rnaseq_normal_sample_name"), 
                        by.y = "sample_name", all.x = TRUE)

#### get the genome hla alleles and number of SNPs ####
genome_hlahd_tables <- output_tables[table_type == "genome_allele_table"]
genome_cohort_alleles <- data.table()
for(line_idx in 1:nrow(genome_hlahd_tables)){
  x <- fread(genome_hlahd_tables[line_idx]$csv_path)
  genome_cohort_alleles <- rbindlist(list(genome_cohort_alleles, x))
}

setnames(genome_cohort_alleles, "num_snps", "n_genome_snps")

# make sure allele 1 is the min (i.e. order them alphabetically)
genome_cohort_alleles[, allele_min := min(c(allele1, allele2)), by = c("patient", "gene")]

genome_cohort_alleles[allele1 == allele_min, new_allele1 := allele1]
genome_cohort_alleles[allele1 == allele_min, new_allele2 := allele2]

genome_cohort_alleles[allele1 != allele_min, new_allele1 := allele2]
genome_cohort_alleles[allele1 != allele_min, new_allele2 := allele1]
genome_cohort_alleles[,c("allele1", "allele2", "allele_min") := NULL]
setnames(genome_cohort_alleles, 
         c("new_allele1", "new_allele2"),
         c("allele1", "allele2"))

#### get the transcriptome hla alleles and number of SNPs ####
transcriptome_hlahd_tables <- output_tables[table_type == "transcriptome_allele_table"]
transcriptome_cohort_alleles <- data.table()
for(line_idx in 1:nrow(transcriptome_hlahd_tables)){
  x <- fread(transcriptome_hlahd_tables[line_idx]$csv_path)
  transcriptome_cohort_alleles <- rbindlist(list(transcriptome_cohort_alleles, x))
}

setnames(transcriptome_cohort_alleles, "num_snps", "n_transcriptome_snps")

# make sure allele 1 is the min (i.e. order them alphabetically)
transcriptome_cohort_alleles[, allele_min := min(c(allele1, allele2)), by = c("patient", "gene")]

transcriptome_cohort_alleles[allele1 == allele_min, new_allele1 := allele1]
transcriptome_cohort_alleles[allele1 == allele_min, new_allele2 := allele2]

transcriptome_cohort_alleles[allele1 != allele_min, new_allele1 := allele2]
transcriptome_cohort_alleles[allele1 != allele_min, new_allele2 := allele1]
transcriptome_cohort_alleles[,c("allele1", "allele2", "allele_min") := NULL]
setnames(transcriptome_cohort_alleles, 
         c("new_allele1", "new_allele2"),
         c("allele1", "allele2"))

# merge the genome and transcriptome together
hlahd_alleles <- merge(genome_cohort_alleles, transcriptome_cohort_alleles,
                       by = c("patient", "gene", "allele1", "allele2", "homozygous"),
                       all = TRUE)

# double check that there weren't different alleles for a given gene and patient
# in the genome or transcritome
if(nrow(hlahd_alleles[,.N,by = c("gene", "patient")][N>1]) > 0){
  stop("Different HLA alleles from the genome and transcriptome allele table.")
}

##### add HLAHD to overview #####
overview_table <- merge(overview_table, 
                        hlahd_alleles, 
                        by = "patient", allow.cartesian=TRUE, all.x = TRUE)

# check each region has an A B and C
n_A <- overview_table[gene == "A", .N, by = "sample_name"]
n_B <- overview_table[gene == "B", .N, by = "sample_name"]
n_C <- overview_table[gene == "C", .N, by = "sample_name"]
setnames(n_A, "N", "A")
setnames(n_B, "N", "B")
setnames(n_C, "N", "C")
overview_table <- merge(overview_table, 
                        n_A, 
                        by = "sample_name", all.x = TRUE)
overview_table <- merge(overview_table, 
                        n_B, 
                        by = "sample_name", all.x = TRUE)
overview_table <- merge(overview_table, 
                        n_C, 
                        by = "sample_name", all.x = TRUE)
overview_table[is.na(A), A:= 0]
overview_table[is.na(B), B:= 0]
overview_table[is.na(C), C:= 0]
if(nrow(overview_table[A != 1 | B != 1 | C != 1])){
  warning("Something has gone wrong with the HLAHD output - there should be 3 genes predictions per patient")
}
overview_table[,c("A", "B", "C") := NULL]

overview_table[, gene := paste0("HLA-", gene)]

##### add in the number of reads before/after bam filtering #####

dna_bam_read_count_tables <- output_tables[table_type == "dna_bam_read_count"]
if(nrow(dna_bam_read_count_tables) > 0){
  cohort_dna_bam_read_count <- data.table()
  for(line_idx in 1:nrow(dna_bam_read_count_tables)){
    x <- fread(dna_bam_read_count_tables[line_idx]$csv_path)
    setnames(x, c("allele", "n_reads_before_filtering", "n_reads_after_filtering"))
    x[,sample_name := gsub("_wxs_novoalign.hla_bam_read_count.csv", "", basename(dna_bam_read_count_tables[line_idx]$csv_path))]
    cohort_dna_bam_read_count <- rbindlist(list(cohort_dna_bam_read_count, x))
  }
  
  # merge in allele 1
  overview_table <- merge(overview_table, cohort_dna_bam_read_count,
                          by.x = c("sample_name", "allele1"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table, 
           c("n_reads_before_filtering", "n_reads_after_filtering"),
           c("allele1_dna_n_reads_before_filtering", "allele1_dna_n_reads_after_filtering"))
  
  # merge in allele 2
  overview_table <- merge(overview_table, cohort_dna_bam_read_count,
                          by.x = c("sample_name", "allele2"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table, 
           c("n_reads_before_filtering", "n_reads_after_filtering"),
           c("allele2_dna_n_reads_before_filtering", "allele2_dna_n_reads_after_filtering"))
  
}else{
  overview_table[,dna_n_reads_before_filtering := NA]
  overview_table[,dna_n_reads_after_filtering := NA]
}

rna_bam_read_count_tables <- output_tables[table_type == "rna_bam_read_count"]
if(nrow(rna_bam_read_count_tables) > 0){
  cohort_rna_bam_read_count <- data.table()
  for(line_idx in 1:nrow(rna_bam_read_count_tables)){
    x <- fread(rna_bam_read_count_tables[line_idx]$csv_path)
    setnames(x, c("allele", "n_reads_before_filtering", "n_reads_after_filtering"))
    x[,sample_name := gsub("_rnaseq_novoalign.hla_bam_read_count.csv", "", basename(rna_bam_read_count_tables[line_idx]$csv_path))]
    cohort_rna_bam_read_count <- rbindlist(list(cohort_rna_bam_read_count, x))
  }
  
  # merge in allele1
  overview_table <- merge(overview_table, cohort_rna_bam_read_count,
                          by.x = c("sample_name", "allele1"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table, 
           c("n_reads_before_filtering", "n_reads_after_filtering"),
           c("allele1_rna_n_reads_before_filtering", "allele1_rna_n_reads_after_filtering"))
  
  # merge in allele2
  overview_table <- merge(overview_table, cohort_rna_bam_read_count,
                          by.x = c("sample_name", "allele2"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table, 
           c("n_reads_before_filtering", "n_reads_after_filtering"),
           c("allele2_rna_n_reads_before_filtering", "allele2_rna_n_reads_after_filtering"))
  
}else{
  overview_table[,rna_n_reads_before_filtering := NA]
  overview_table[,rna_n_reads_after_filtering := NA]
}

#### add in the dna analysis ####
dna_tables <- output_tables[table_type == "dna_analysis"]
if(nrow(dna_tables) > 0){
  cohort_dna <- data.table()
  for(line_idx in 1:nrow(dna_tables)){
    x <- fread(dna_tables[line_idx]$csv_path)
    
    allele1_dt <- x[,c("sample_name", "allele1", 
                       "cn1_binned", "cn1_binned_lower", "cn1_binned_upper",
                       "cn_n_bins", "cn_n_snps",
                       "allele1_exp_dp",
                       "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test",
                       "logr_aib_n_snps")]
    allele2_dt <- x[,c("sample_name", "allele2", 
                       "cn2_binned", "cn2_binned_lower", "cn2_binned_upper",
                       "cn_n_bins", "cn_n_snps",
                       "allele2_exp_dp",
                       "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test",
                       "logr_aib_n_snps")]
    
    allele_dt <- rbindlist(list(allele1_dt, allele2_dt), use.names = FALSE)
    
    cohort_dna <- rbindlist(list(cohort_dna, allele_dt))
  }
  
  # merge 
  setnames(cohort_dna, 
           old = c("allele1", "cn1_binned", "cn1_binned_lower", "cn1_binned_upper", "allele1_exp_dp"), 
           new = c("allele", "cn_binned", "cn_binned_lower", "cn_binned_upper", "allele_expected_depth"))
  
  # merge allele 1
  overview_table <- merge(overview_table, 
                          cohort_dna, 
                          by.x = c("sample_name", "allele1"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("cn_binned", "cn_binned_lower", "cn_binned_upper", "allele_expected_depth"),
           c("cn1_binned", "cn1_binned_lower", "cn1_binned_upper", "allele1_expected_depth"))
  
  cohort_dna[,c("cn_n_bins", "cn_n_snps", "logr_aib_paired_t_test", 
                "logr_aib_paired_wilcoxon_test", "logr_aib_n_snps") := NULL]
  
  # merge allele 2
  overview_table <- merge(overview_table, 
                          cohort_dna, 
                          by.x = c("sample_name", "allele2"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("cn_binned", "cn_binned_lower", "cn_binned_upper", "allele_expected_depth"),
           c("cn2_binned", "cn2_binned_lower", "cn2_binned_upper", "allele2_expected_depth"))
}else{
  overview_table[,cn1_binned := NA]
  overview_table[,cn1_binned_lower := NA]
  overview_table[,cn1_binned_upper := NA]
  overview_table[,allele1_expected_depth := NA]
  overview_table[,cn_n_bins := NA]
  overview_table[,cn_n_snps := NA]
  overview_table[,logr_aib_paired_t_test := NA]
  overview_table[,logr_aib_paired_wilcoxon_test := NA]
  overview_table[,logr_aib_n_snps := NA]
  
  overview_table[,cn2_binned := NA]
  overview_table[,cn2_binned_lower := NA]
  overview_table[,cn2_binned_upper := NA]
  overview_table[,allele2_expected_depth := NA]
}

#### Add rpkm ####
rpkm_tables <- output_tables[table_type == "rpkm"]
if(nrow(rpkm_tables) > 0){
  cohort_rpkm <- data.table()
  for(line_idx in 1:nrow(rpkm_tables)){
    x <- fread(rpkm_tables[line_idx]$csv_path)
    
    allele1_dt <- x[,c( "sample_name", "gene_rpkm", "allele1", 
                        "allele1_rpkm", 
                        "frac_mapping_multi_gene", "frac_mapping_uniquely")]
    allele2_dt <- x[,c("sample_name", "gene_rpkm", "allele2", 
                       "allele2_rpkm", 
                       "frac_mapping_multi_gene", "frac_mapping_uniquely")]
    
    allele_dt <- rbindlist(list(allele1_dt, allele2_dt), use.names = FALSE)
    setnames(allele_dt, 
             c("allele1", "allele1_rpkm"),
             c("allele", "allele_rpkm"))
    # remove empty allels due to homozygous samples
    allele_dt <- allele_dt[allele != ""]
    cohort_rpkm <- rbindlist(list(cohort_rpkm, allele_dt))
  }
  
  # merge allele1
  overview_table <- merge(overview_table, 
                          cohort_rpkm, 
                          by.x = c("sample_name", "allele1"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("allele_rpkm"),
           c("allele1_rpkm"))
  cohort_rpkm[,c("gene_rpkm", "frac_mapping_multi_gene", "frac_mapping_uniquely") := NULL]
  
  # merge allele2
  overview_table <- merge(overview_table, 
                          cohort_rpkm, 
                          by.x = c("sample_name", "allele2"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("allele_rpkm"),
           c("allele2_rpkm"))
  
  overview_table[homozygous == TRUE, allele2_rpkm := NA]
  
}else{
  overview_table[,gene_rpkm := NA]
  overview_table[,allele1_rpkm := NA]
  overview_table[,allele2_rpkm := NA]
  overview_table[,frac_mapping_multi_gene := NA]
  overview_table[,frac_mapping_uniquely := NA]
}

#### Add in rna aib ####
rna_aib_tables <- output_tables[table_type == "rna_aib"]
if(nrow(rna_aib_tables) > 0){
  cohort_rna_aib <- data.table()
  for(line_idx in 1:nrow(rna_aib_tables)){
    x <- fread(rna_aib_tables[line_idx]$csv_path)
    cohort_rna_aib <- rbindlist(list(cohort_rna_aib, x))
  }
  
  cohort_rna_aib[,c("allele1", "allele2", "rna_aib_allele1_n_snps_with_coverage", "rna_aib_allele2_n_snps_with_coverage") := NULL]
  
  # merge
  overview_table <- merge(overview_table,
                          cohort_rna_aib,
                          by = c("sample_name", "gene"),
                          all.x = TRUE)
  
}else{
  overview_table[,rna_aib_paired_t_test := NA]
  overview_table[,rna_aib_paired_wilcoxon_test := NA]
}

#### add in the tumour normal p value ####
rna_tumour_normal_tables <- output_tables[table_type == "repression"]
if(nrow(rna_tumour_normal_tables) > 0){
  cohort_tumour_normal <- data.table()
  for(line_idx in 1:nrow(rna_tumour_normal_tables)){
    x <- fread(rna_tumour_normal_tables[line_idx]$csv_path)
    cohort_tumour_normal <- rbindlist(list(cohort_tumour_normal,x))
  }
  cohort_tumour_normal[,c("normal_sample_name", "repression_normal_n_snps_with_coverage",
                          "repression_tumour_n_snps_with_coverage") := NULL]
  
  # merge allele1
  overview_table <- merge(overview_table, 
                          cohort_tumour_normal, 
                          by.x = c("sample_name", "allele1"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("repression_paired_t_test", "repression_paired_wilcoxon_test",
             "repression_median_tumour_dp", "repression_median_normal_dp"),
           c("allele1_repression_paired_t_test", "allele1_repression_paired_wilcoxon_test",
             "allele1_repression_median_tumour_dp", "allele1_repression_median_normal_dp"))
  
  # merge allele2
  overview_table <- merge(overview_table, 
                          cohort_tumour_normal, 
                          by.x = c("sample_name", "allele2"),
                          by.y = c("sample_name", "allele"),
                          all.x = TRUE)
  setnames(overview_table,
           c("repression_paired_t_test", "repression_paired_wilcoxon_test",
             "repression_median_tumour_dp", "repression_median_normal_dp"),
           c("allele2_repression_paired_t_test", "allele2_repression_paired_wilcoxon_test",
             "allele2_repression_median_tumour_dp", "allele2_repression_median_normal_dp"))
  
}else{
  overview_table[,allele1_repression_paired_t_test := NA]
  overview_table[,allele1_repression_paired_wilcoxon_test := NA]
  overview_table[,allele1_repression_median_tumour_dp := NA]
  overview_table[,allele1_repression_median_normal_dp := NA]
  overview_table[,allele2_repression_paired_t_test := NA]
  overview_table[,allele2_repression_paired_wilcoxon_test := NA]
  overview_table[,allele2_repression_median_tumour_dp := NA]
  overview_table[,allele2_repression_median_normal_dp := NA]
}

###### fail samples ######

# homozygous
overview_table[homozygous == TRUE, fail_homozygous := TRUE]
overview_table[homozygous == FALSE, fail_homozygous := FALSE]

# fail if too few dna snps
overview_table[sample_has_wes == TRUE & 
                 (cn_n_snps < min_n_snps | logr_aib_n_snps < min_n_snps | fail_homozygous == TRUE), fail_n_wes_snps := TRUE]
overview_table[sample_has_wes == TRUE & 
                 (cn_n_snps >= min_n_snps & logr_aib_n_snps >= min_n_snps), fail_n_wes_snps := FALSE]

# overview_table[sample_has_wes == TRUE,.N,by = fail_n_wes_snps]

# fail if low expected depth
overview_table[sample_has_wes == TRUE & 
                 (allele1_expected_depth < min_expected_depth | allele2_expected_depth < min_expected_depth | fail_homozygous == TRUE | cn_n_snps == 0), 
               fail_expected_depth := TRUE]
overview_table[sample_has_wes == TRUE & (allele1_expected_depth >= min_expected_depth & allele2_expected_depth >= min_expected_depth), 
               fail_expected_depth := FALSE]

overview_table[sample_has_wes == TRUE,.N,by = fail_expected_depth]

# cn range
overview_table[sample_has_wes == TRUE, allele1_cn_range := cn1_binned_upper - cn1_binned_lower]
overview_table[sample_has_wes == TRUE, allele2_cn_range := cn2_binned_upper - cn2_binned_lower]

overview_table[sample_has_wes == TRUE & (allele1_cn_range > max_cn_range | allele2_cn_range > max_cn_range),
               fail_cn_range := TRUE]
overview_table[sample_has_wes == TRUE & (is.na(allele1_cn_range) | is.na(allele2_cn_range)), 
               fail_cn_range := TRUE]
overview_table[sample_has_wes == TRUE & allele1_cn_range <= max_cn_range & allele2_cn_range <= max_cn_range,
               fail_cn_range := FALSE]

overview_table[sample_has_wes == TRUE,.N,by = fail_cn_range]

# fail if too few rna aib snps
overview_table[sample_has_rna == TRUE & n_transcriptome_snps < min_n_snps, fail_rna_n_snps := TRUE]
overview_table[sample_has_rna == TRUE & n_transcriptome_snps >= min_n_snps, fail_rna_n_snps := FALSE]

# overview_table[,.N,by = c("fail_rna_n_snps", "sample_has_rna")]

# fail if too many reads map to multi gene
overview_table[sample_has_rna == TRUE & frac_mapping_multi_gene > max_frac_mapping_multi_gene , 
               fail_multi_gene_mapping := TRUE]
overview_table[sample_has_rna == TRUE & frac_mapping_multi_gene <= max_frac_mapping_multi_gene , 
               fail_multi_gene_mapping := FALSE]

# overview_table[,.N,by = c("fail_multi_gene_mapping", "sample_has_rna")]

# fail if too many reads map to both alleles
overview_table[sample_has_rna == TRUE & frac_mapping_uniquely > min_frac_mapping_uniquely , 
               fail_unique_allele_mapping := FALSE]
overview_table[sample_has_rna == TRUE & frac_mapping_uniquely <= min_frac_mapping_uniquely , 
               fail_unique_allele_mapping := TRUE]
overview_table[,.N,by = c("fail_unique_allele_mapping", "sample_has_rna")]

# pass/fail the wes
overview_table[sample_has_wes == TRUE &
                 (fail_homozygous == TRUE | 
                    fail_n_wes_snps == TRUE |
                    fail_expected_depth == TRUE |
                    fail_cn_range == TRUE), wes_fail := TRUE]
overview_table[sample_has_wes == TRUE &
                 (fail_homozygous == FALSE & 
                    fail_n_wes_snps == FALSE &
                    fail_expected_depth == FALSE & 
                    fail_cn_range == FALSE), wes_fail := FALSE]

# pass/fail the rnaseq
overview_table[sample_has_rna == TRUE &
                 (fail_homozygous == TRUE | 
                    fail_rna_n_snps == TRUE |
                    fail_multi_gene_mapping == TRUE | 
                    fail_unique_allele_mapping == TRUE), rna_fail := TRUE]
overview_table[sample_has_rna == TRUE &
                 (fail_homozygous == FALSE & 
                    fail_rna_n_snps == FALSE &
                    fail_multi_gene_mapping == FALSE & 
                    fail_unique_allele_mapping == FALSE), rna_fail := FALSE]

if(nrow(overview_table[rnaseq_normal_sample_name != ""]) > 0 & "rna_fail" %in% colnames(overview_table)){
  rna_normal_fail_dt <- overview_table[sample_type == "normal" & sample_has_rna == TRUE,
                                       c("sample_name", "gene", "rna_fail")]
  setnames(rna_normal_fail_dt, "rna_fail", "rna_normal_fail")
  overview_table <- merge(overview_table, rna_normal_fail_dt, 
                          by.x = c("rnaseq_normal_sample_name", "gene"),
                          by.y = c("sample_name", "gene"), all.x = TRUE)
  
  overview_table[sample_has_rna == TRUE & rnaseq_normal_sample_name != "" &
                   (rna_normal_fail == TRUE | rna_fail == TRUE) , rna_matched_normal_fail := TRUE]
  overview_table[sample_has_rna == TRUE & rnaseq_normal_sample_name != "" &
                   (rna_normal_fail == FALSE & rna_fail == FALSE) , rna_matched_normal_fail := FALSE]
}else{
  overview_table[,rna_normal_fail := NA]
  overview_table[,rna_repression_fail := NA]
}

# add in the matched normal rna rpkm
normal_rpkm <- overview_table[rna_fail == FALSE & sample_type == "normal", 
                     c("sample_name", "gene",
                       "allele1_rpkm", "allele2_rpkm")]
setnames(normal_rpkm, 
         c("allele1_rpkm", "allele2_rpkm"), 
         c("allele1_matched_normal_rpkm", "allele2_matched_normal_rpkm"))
overview_table <- merge(overview_table, normal_rpkm, 
                        by.x = c("rnaseq_normal_sample_name", "gene"),
                        by.y = c("sample_name", "gene"), all.x = TRUE)

# update the allele name
overview_table[,allele1 := toupper(allele1)]
overview_table[,allele1 := gsub("HLA_", "HLA-", allele1)]
overview_table[,allele1 := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele1)]
overview_table[,allele1 := gsub("_", ":", allele1)]

overview_table[,allele2 := toupper(allele2)]
overview_table[,allele2 := gsub("HLA_", "HLA-", allele2)]
overview_table[,allele2 := gsub("(^.*?)_(.*)", "\\1\\*\\2", allele2)]
overview_table[,allele2 := gsub("_", ":", allele2)]

# add in the short allele
overview_table[, short_allele1 := gsub("(HLA-.*?:.*?):.*", "\\1", allele1)]
overview_table[, short_allele2 := gsub("(HLA-.*?:.*?):.*", "\\1", allele2)]

# add in dna aib
overview_table[wes_fail == FALSE & logr_aib_paired_wilcoxon_test < 0.01, dna_aib := TRUE]
overview_table[wes_fail == FALSE & logr_aib_paired_wilcoxon_test >= 0.01, dna_aib := FALSE]
# overview_table[,.N,by = c("wes_fail", "dna_aib", "sample_type")]

# add in loh
overview_table[wes_fail == FALSE & dna_aib == TRUE & cn1_binned < 0.5, allele1_loh := TRUE]
overview_table[wes_fail == FALSE & (dna_aib == FALSE | cn1_binned >= 0.5), allele1_loh := FALSE]

overview_table[wes_fail == FALSE & dna_aib == TRUE & cn2_binned < 0.5, allele2_loh := TRUE]
overview_table[wes_fail == FALSE & (dna_aib == FALSE | cn2_binned >= 0.5), allele2_loh := FALSE]

overview_table[allele1_loh == TRUE | allele2_loh == TRUE, loh := TRUE]
overview_table[allele1_loh == FALSE & allele2_loh == FALSE, loh := FALSE]

overview_table[wes_fail == FALSE & cn1_binned < cn2_binned, minor_allele := allele1]
overview_table[wes_fail == FALSE & cn2_binned < cn1_binned, minor_allele := allele2]

# overview_table[,.N,by = c("wes_fail", "allele1_loh", "sample_type")]
# overview_table[,.N,by = c("wes_fail", "allele2_loh", "sample_type")]
# overview_table[,.N,by = c("wes_fail", "loh", "sample_type")]

# add in rna aib
overview_table[rna_fail == FALSE & rna_aib_paired_wilcoxon_test < 0.01, rna_aib := TRUE]
overview_table[rna_fail == FALSE & rna_aib_paired_wilcoxon_test >= 0.01, rna_aib := FALSE]

# add in allele1 repression
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 allele1_repression_paired_wilcoxon_test < 0.01 &
                 allele1_repression_median_tumour_dp < allele1_repression_median_normal_dp,
               allele1_repressed := TRUE]
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 (allele1_repression_paired_wilcoxon_test >= 0.01 |
                    allele1_repression_median_tumour_dp > allele1_repression_median_normal_dp),
               allele1_repressed := FALSE]

overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 allele1_repression_paired_wilcoxon_test < 0.01 &
                 allele1_repression_median_tumour_dp > allele1_repression_median_normal_dp,
               allele1_over_expressed := TRUE]
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 (allele1_repression_paired_wilcoxon_test >= 0.01 |
                    allele1_repression_median_tumour_dp < allele1_repression_median_normal_dp),
               allele1_over_expressed := FALSE]

# add in allele2 repression
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 allele2_repression_paired_wilcoxon_test < 0.01 &
                 allele2_repression_median_tumour_dp < allele2_repression_median_normal_dp,
               allele2_repressed := TRUE]
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 (allele2_repression_paired_wilcoxon_test >= 0.01 |
                    allele2_repression_median_tumour_dp > allele2_repression_median_normal_dp),
               allele2_repressed := FALSE]

overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 allele2_repression_paired_wilcoxon_test < 0.01 &
                 allele2_repression_median_tumour_dp > allele2_repression_median_normal_dp,
               allele2_over_expressed := TRUE]
overview_table[rna_fail == FALSE & rna_normal_fail == FALSE & 
                 (allele2_repression_paired_wilcoxon_test >= 0.01 |
                    allele2_repression_median_tumour_dp < allele2_repression_median_normal_dp),
               allele2_over_expressed := FALSE]


# overview_table[,.N,by = c("rna_fail", "rna_normal_fail", "allele2_repressed")]
# overview_table[,.N,by = c("rna_fail", "rna_normal_fail", "allele1_repressed")]
# overview_table[,.N,by = c("rna_fail", "rna_normal_fail", "allele1_over_expressed")]
# overview_table[,.N,by = c("rna_fail", "rna_normal_fail", "allele2_over_expressed")]

overview_table[, wes_min_snp_depth := dna_snp_min_depth]

# order the columns
col_order <- c("patient", "sample_name", "sample_type", "gene", "allele1", "allele2",
               "short_allele1", "short_allele2",
               "homozygous", "purity", "ploidy",
               "sample_has_wes", "sample_has_rna", "rna_fail", "wes_fail",
               "wes_germline_id", "rnaseq_normal_sample_name",
               "n_genome_snps", "n_transcriptome_snps", 
               "allele1_dna_n_reads_before_filtering", "allele1_dna_n_reads_after_filtering", "allele2_dna_n_reads_before_filtering", "allele2_dna_n_reads_after_filtering",
               "allele1_rna_n_reads_before_filtering", "allele1_rna_n_reads_after_filtering", "allele2_rna_n_reads_before_filtering", "allele2_rna_n_reads_after_filtering", 
               "cn1_binned", "cn1_binned_lower", "cn1_binned_upper", "allele1_cn_range",
               "cn2_binned", "cn2_binned_lower", "cn2_binned_upper", "allele2_cn_range",
               "cn_n_bins", "cn_n_snps", "wes_min_snp_depth",
               "allele1_expected_depth", "allele2_expected_depth", 
               "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test", "logr_aib_n_snps", 
               "dna_aib", "loh", "minor_allele",
               "gene_rpkm", "allele1_rpkm", "allele2_rpkm", 
               "rna_aib_paired_t_test", "rna_aib_paired_wilcoxon_test",
               "rna_aib",
               "frac_mapping_multi_gene", "frac_mapping_uniquely", 
               "allele1_repression_paired_t_test", "allele1_repression_paired_wilcoxon_test", 
               "allele1_repression_median_tumour_dp", "allele1_repression_median_normal_dp", 
               "allele2_repression_paired_t_test", "allele2_repression_paired_wilcoxon_test", 
               "allele2_repression_median_tumour_dp", "allele2_repression_median_normal_dp", 
               "allele1_repressed", "allele1_over_expressed", "allele2_repressed", "allele2_over_expressed",
               "fail_homozygous", "fail_n_wes_snps", "fail_expected_depth", "fail_cn_range",
               "fail_rna_n_snps", "fail_multi_gene_mapping", "fail_unique_allele_mapping")

setcolorder(overview_table, col_order)

# make allele level table
allele1_dt <- overview_table[,c("patient", "sample_name", "sample_type", "gene", "allele1", "short_allele1",
                                "homozygous", "purity", "ploidy",
                                "sample_has_wes", "sample_has_rna", "rna_fail", "wes_fail", "rna_normal_fail",
                                "wes_germline_id", "rnaseq_normal_sample_name",
                                "n_genome_snps", "n_transcriptome_snps", 
                                "allele1_dna_n_reads_before_filtering", "allele1_dna_n_reads_after_filtering", 
                                "allele1_rna_n_reads_before_filtering", "allele1_rna_n_reads_after_filtering", 
                                "cn1_binned", "cn1_binned_lower", "cn1_binned_upper", "allele1_cn_range",
                                "cn_n_bins", "cn_n_snps", "wes_min_snp_depth",
                                "allele1_expected_depth", 
                                "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test", "logr_aib_n_snps", 
                                "dna_aib", "allele1_loh", "minor_allele",
                                "gene_rpkm", "allele1_rpkm", "allele1_matched_normal_rpkm",
                                "rna_aib_paired_t_test", "rna_aib_paired_wilcoxon_test",
                                "rna_aib",
                                "frac_mapping_multi_gene", "frac_mapping_uniquely", 
                                "allele1_repression_paired_t_test", "allele1_repression_paired_wilcoxon_test", 
                                "allele1_repression_median_tumour_dp", "allele1_repression_median_normal_dp", 
                                "allele1_repressed", "allele1_over_expressed", 
                                "fail_homozygous", "fail_n_wes_snps", "fail_expected_depth", "fail_cn_range",
                                "fail_rna_n_snps", "fail_multi_gene_mapping", "fail_unique_allele_mapping")]

setnames(allele1_dt,
         old = c("allele1", "short_allele1",
                 "allele1_dna_n_reads_before_filtering", "allele1_dna_n_reads_after_filtering", 
                 "allele1_rna_n_reads_before_filtering", "allele1_rna_n_reads_after_filtering", 
                 "cn1_binned", "cn1_binned_lower", "cn1_binned_upper", "allele1_cn_range",
                 "allele1_expected_depth", "allele1_loh",
                 "allele1_rpkm", "allele1_matched_normal_rpkm", 
                 "allele1_repression_paired_t_test", "allele1_repression_paired_wilcoxon_test", 
                 "allele1_repression_median_tumour_dp", "allele1_repression_median_normal_dp", 
                 "allele1_repressed", "allele1_over_expressed"),
         new = c("allele", "short_allele",
                 "dna_n_reads_before_filtering", "dna_n_reads_after_filtering", 
                 "rna_n_reads_before_filtering", "rna_n_reads_after_filtering", 
                 "cn_binned", "cn_binned_lower", "cn_binned_upper", "cn_range",
                 "expected_depth", "loh",
                 "allele_rpkm", "allele_matched_normal_rpkm",
                 "repression_paired_t_test", "repression_paired_wilcoxon_test", 
                 "repression_median_tumour_dp", "repression_median_normal_dp", 
                 "repressed", "over_expressed"))

allele2_dt <- overview_table[homozygous == FALSE,
                             c("patient", "sample_name", "sample_type", "gene", "allele2", "short_allele2",
                               "homozygous", "purity", "ploidy",
                               "sample_has_wes", "sample_has_rna", "rna_fail", "wes_fail", "rna_normal_fail",
                               "wes_germline_id", "rnaseq_normal_sample_name",
                               "n_genome_snps", "n_transcriptome_snps", 
                               "allele2_dna_n_reads_before_filtering", "allele2_dna_n_reads_after_filtering", 
                               "allele2_rna_n_reads_before_filtering", "allele2_rna_n_reads_after_filtering", 
                               "cn2_binned", "cn2_binned_lower", "cn2_binned_upper", "allele2_cn_range",
                               "cn_n_bins", "cn_n_snps", "wes_min_snp_depth",
                               "allele2_expected_depth", 
                               "logr_aib_paired_t_test", "logr_aib_paired_wilcoxon_test", "logr_aib_n_snps", 
                               "dna_aib", "allele2_loh", "minor_allele",
                               "gene_rpkm", "allele2_rpkm", "allele2_matched_normal_rpkm",
                               "rna_aib_paired_t_test", "rna_aib_paired_wilcoxon_test",
                               "rna_aib",
                               "frac_mapping_multi_gene", "frac_mapping_uniquely", 
                               "allele2_repression_paired_t_test", "allele2_repression_paired_wilcoxon_test", 
                               "allele2_repression_median_tumour_dp", "allele2_repression_median_normal_dp", 
                               "allele2_repressed", "allele2_over_expressed", 
                               "fail_homozygous", "fail_n_wes_snps", "fail_expected_depth", "fail_cn_range",
                               "fail_rna_n_snps", "fail_multi_gene_mapping", "fail_unique_allele_mapping")]

setnames(allele2_dt,
         old = c("allele2", "short_allele2",
                 "allele2_dna_n_reads_before_filtering", "allele2_dna_n_reads_after_filtering", 
                 "allele2_rna_n_reads_before_filtering", "allele2_rna_n_reads_after_filtering", 
                 "cn2_binned", "cn2_binned_lower", "cn2_binned_upper", "allele2_cn_range",
                 "allele2_expected_depth", "allele2_loh",
                 "allele2_rpkm", "allele2_matched_normal_rpkm",
                 "allele2_repression_paired_t_test", "allele2_repression_paired_wilcoxon_test", 
                 "allele2_repression_median_tumour_dp", "allele2_repression_median_normal_dp", 
                 "allele2_repressed", "allele2_over_expressed"),
         new = c("allele", "short_allele",
                 "dna_n_reads_before_filtering", "dna_n_reads_after_filtering", 
                 "rna_n_reads_before_filtering", "rna_n_reads_after_filtering", 
                 "cn_binned", "cn_binned_lower", "cn_binned_upper", "cn_range",
                 "expected_depth", "loh",
                 "allele_rpkm", "allele_matched_normal_rpkm",
                 "repression_paired_t_test", "repression_paired_wilcoxon_test", 
                 "repression_median_tumour_dp", "repression_median_normal_dp", 
                 "repressed", "over_expressed"))

allele_table <- rbindlist(list(allele1_dt, allele2_dt))

# double check that neither table has extra rows
if(nrow(allele_table[,.N,by = c("sample_name", "allele")][N>1]) > 0){
  stop("allele table has more than one row per sample/allele")
}

if(nrow(overview_table[,.N,by = c("sample_name", "gene")][N>1]) > 0){
  stop("allele table has more than one row per sample/allele")
}

fwrite(overview_table, "cohort_mhc_hammer_gene_table.csv")
fwrite(allele_table, "cohort_mhc_hammer_allele_table.csv")

