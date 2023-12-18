suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--vep_tables',
                    help='Path to VEP txt files.',
                    nargs="+",
                    required=TRUE)

parser$add_argument('--wxs_tumour_bam_files',
                    help='Path to WXS tumour bams.',
                    nargs="+",
                    required=TRUE)

parser$add_argument('--wxs_gl_bam_files',
                    help='Sample names of germlines.',
                    nargs="+",
                    required=TRUE)

parser$add_argument('--mutation_save_path',
                    help='Path to save mutation file as',
                    required=TRUE)

parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

parser$add_argument('--inventory_path',  nargs=1,
                    help='Path to cohort inventory',
                    required=TRUE)


args <- parser$parse_args()
vep_tables <- args$vep_tables
wxs_tumour_bam_files <- args$wxs_tumour_bam_files
wxs_gl_bam_files <- args$wxs_gl_bam_files
mutation_save_path <- args$mutation_save_path
scripts_dir <- args$scripts_dir
inventory <- fread(args$inventory_path)

source(paste0(scripts_dir, "/mutation_table_function.R"))

cat("vep_tables=", vep_tables, "\n")
cat("wxs_tumour_bam_files=", wxs_tumour_bam_files, "\n")
cat("wxs_gl_bam_files=", wxs_gl_bam_files, "\n")
cat("mutation_save_path=", mutation_save_path, "\n")
cat("scripts_dir=", scripts_dir, "\n")
cat("inventory_path=", args$inventory_path, "\n")

tumour_samples <- unique(gsub("_wxs_novoalign.*", "", wxs_tumour_bam_files))
germline_samples <- unique(gsub("_wxs_novoalign.*", "", wxs_gl_bam_files))

# first get a unique set of all mutations
tumour_mutations <- data.table()
for(vep_tab in vep_tables){
  
  region_mutations <- fread(vep_tab)
  
  tumour_dp_col <- names(region_mutations)[names(region_mutations) %in% paste0(tumour_samples, ".DP")]
  tumour_ad_col <- names(region_mutations)[names(region_mutations) %in% paste0(tumour_samples, ".AD")]
  
  germline_dp_col <- names(region_mutations)[names(region_mutations) %in% paste0(germline_samples, ".DP")]
  germline_ad_col <- names(region_mutations)[names(region_mutations) %in% paste0(germline_samples, ".AD")]
  
  if(length(tumour_dp_col) != 1 |
     length(tumour_ad_col) != 1 |
     length(germline_dp_col) != 1 |
     length(germline_ad_col) != 1 ){
    stop("There should be a tumour and germline column")
  }
  
  setnames(region_mutations, 
           old = c(tumour_dp_col, tumour_ad_col, germline_dp_col, germline_ad_col),
           new = c("tumour_dp", "tumour_ad", "germline_dp", "germline_ad"))
  
  region_mutations[,tumour_sample_name := gsub(".DP$", "", tumour_dp_col)]
  
  tumour_mutations <- rbindlist(list(tumour_mutations, region_mutations), use.names = TRUE)
  
}

tumour_mutations[,tumour_ref_dp := gsub(",.*", "", tumour_ad)]
tumour_mutations[,tumour_alt_dp := gsub(".*,", "", tumour_ad)]

tumour_mutations[,germline_ref_dp := gsub(",.*", "", germline_ad)]
tumour_mutations[,germline_alt_dp := gsub(".*,", "", germline_ad)]

# split the VEP consequence line in separate columns
for(line_idx in 1:nrow(tumour_mutations)){
  vep_line <- tumour_mutations[line_idx]$CSQ
  vep_line_split <- strsplit(x = vep_line,  split =  "|", fixed = TRUE)[[1]]
  tumour_mutations[line_idx, vep_impact := vep_line_split[3]]
  tumour_mutations[line_idx, vep_feature_type := vep_line_split[4]]
  tumour_mutations[line_idx, vep_feature := vep_line_split[5]]
  tumour_mutations[line_idx, vep_exon := vep_line_split[6]]
  tumour_mutations[line_idx, vep_intron := vep_line_split[7]]
  tumour_mutations[line_idx, vep_cdna_position := vep_line_split[8]]
  tumour_mutations[line_idx, vep_cds_position := vep_line_split[9]]
  tumour_mutations[line_idx, vep_protein_position := vep_line_split[10]]
  tumour_mutations[line_idx, vep_amino_acids := vep_line_split[11]]
  tumour_mutations[line_idx, vep_codons := vep_line_split[12]]
  tumour_mutations[line_idx, vep_existing_variation := vep_line_split[13]]
  tumour_mutations[line_idx, vep_existing_distance := vep_line_split[14]]
  tumour_mutations[line_idx, vep_existing_strand := vep_line_split[15]]
  tumour_mutations[line_idx, vep_existing_flag := vep_line_split[16]]
  
  # vep consequence we want to split and order alphabetically
  vep_consequence_line <- vep_line_split[2]
  vep_consequence_split <- strsplit(vep_consequence_line, "&")[[1]]
  vep_consequence_split <- sort(vep_consequence_split)
  new_vep_consequence <- paste0(vep_consequence_split, collapse = "&")
  tumour_mutations[line_idx, vep_consequence := new_vep_consequence]
}

# add in mutation type
tumour_mutations[,nchar_ref := nchar(REF)]
tumour_mutations[,nchar_alt := nchar(ALT)]
tumour_mutations[nchar_ref == 1 & nchar_alt == 1, mut_type := "SNV"]
tumour_mutations[nchar_ref > nchar_alt, mut_type := "DEL"]
tumour_mutations[nchar_ref < nchar_alt, mut_type := "INS"]
tumour_mutations[nchar_ref == nchar_alt & nchar_alt > 1, mut_type := "MNV"]

# add in the start and end dna position
tumour_mutations[,start_dna := POS]
tumour_mutations[mut_type == "MNV", end_dna := POS + nchar_alt - 1]
tumour_mutations[mut_type == "DEL", end_dna := POS + nchar_ref - 1]

tumour_mutations[,c("nchar_ref", "nchar_alt") := NULL]

# loop over the tumour samples and get the ref and var counts on the unique set of mutations

unique_mutations <- unique(tumour_mutations[,c("CHROM", "POS", "REF", "ALT",
                                               "mut_type", "start_dna", "end_dna")])
tumour_bam_read_counts <- data.table()
for(tumour_sample in tumour_samples){
  
  region_bam_read_counts <- copy(unique_mutations)
  region_bam_read_counts[,tumour_sample_name := tumour_sample]
  region_bam_read_counts[,bam_path := paste0(tumour_sample, "_wxs_novoalign.",
                                             CHROM, ".sorted.filtered.bam")]
  
  for(line_idx in 1:nrow(region_bam_read_counts)){
    
    mutation_type <- region_bam_read_counts[line_idx]$mut_type
    chr <- region_bam_read_counts[line_idx]$CHROM
    start <- region_bam_read_counts[line_idx]$start_dna
    stop <- region_bam_read_counts[line_idx]$end_dna
    ref <- region_bam_read_counts[line_idx]$REF
    alt <- region_bam_read_counts[line_idx]$ALT
    
    bam_read_count_output <- get_bam_read_count(bam_path = region_bam_read_counts[line_idx]$bam_path, 
                                                chr, start, stop, ref, alt, mutation_type)
    
    region_bam_read_counts[line_idx, tumour_ref_count := bam_read_count_output$bam_ref_count]
    region_bam_read_counts[line_idx, tumour_alt_count := bam_read_count_output$bam_alt_count]
    region_bam_read_counts[line_idx, tumour_N_count := bam_read_count_output$bam_N_count]
    
  }
  
  tumour_bam_read_counts <- rbindlist(list(tumour_bam_read_counts, region_bam_read_counts))
}

germline_bam_read_counts <- data.table()
for(germline_sample in germline_samples){
  
  region_bam_read_counts <- copy(unique_mutations)
  region_bam_read_counts[,germline_sample_name := germline_sample]
  region_bam_read_counts[,bam_path := paste0(germline_sample, "_wxs_novoalign.",
                                             CHROM, ".sorted.filtered.bam")]
  
  for(line_idx in 1:nrow(region_bam_read_counts)){
    
    mutation_type <- region_bam_read_counts[line_idx]$mut_type
    chr <- region_bam_read_counts[line_idx]$CHROM
    start <- region_bam_read_counts[line_idx]$start_dna
    stop <- region_bam_read_counts[line_idx]$end_dna
    ref <- region_bam_read_counts[line_idx]$REF
    alt <- region_bam_read_counts[line_idx]$ALT
    
    bam_read_count_output <- get_bam_read_count(bam_path = region_bam_read_counts[line_idx]$bam_path, 
                                                chr, start, stop, ref, alt, mutation_type)
    
    region_bam_read_counts[line_idx, germline_ref_count := bam_read_count_output$bam_ref_count]
    region_bam_read_counts[line_idx, germline_alt_count := bam_read_count_output$bam_alt_count]
    region_bam_read_counts[line_idx, germline_N_count := bam_read_count_output$bam_N_count]
    
  }
  
  germline_bam_read_counts <- rbindlist(list(germline_bam_read_counts, region_bam_read_counts))
}

# add the matched normal names to the tumour table
tumour_normal_mapping <- unique(inventory[sample_type == "tumour" & 
                                            sequencing_type == "wxs" &
                                            sample_name %in% tumour_bam_read_counts$tumour_sample_name,
                                          c("sample_name", "normal_sample_name")])
setnames(tumour_normal_mapping, "normal_sample_name", "germline_sample_name")

tumour_bam_read_counts <- merge(tumour_bam_read_counts, tumour_normal_mapping,
                                by.x = "tumour_sample_name", by.y = "sample_name",
                                all.x = TRUE)

if(nrow(tumour_bam_read_counts[is.na(germline_sample_name)]) > 0){
  stop("Why are some tumours missing matched germlines?")
}

# add the normal bam read counts to the tumour table
germline_bam_read_counts <- germline_bam_read_counts[,c("CHROM", "POS", "REF", "ALT", 
                                                        "germline_sample_name", 
                                                        "germline_ref_count",
                                                        "germline_alt_count",
                                                        "germline_N_count")]
tumour_bam_read_counts <- merge(tumour_bam_read_counts, germline_bam_read_counts,
                                by = c("CHROM", "POS", "REF", "ALT", "germline_sample_name"))

# add back the VEP CSQs ect
tumour_mutations[,c("mut_type", "start_dna", "end_dna") := NULL]
tumour_bam_read_counts[,c("start_dna", "end_dna") := NULL]
tumour_bam_read_counts <- merge(tumour_bam_read_counts, tumour_mutations,
                                by = c("CHROM", "POS", "REF", "ALT", "tumour_sample_name"),
                                all.x = TRUE)

# add in if it was called by mutect
tumour_bam_read_counts[is.na(FILTER), called_by_mutect := FALSE]
tumour_bam_read_counts[!is.na(FILTER), called_by_mutect := TRUE]

setnames(tumour_bam_read_counts, 
         c("CHROM", "POS", "REF", "ALT", "FILTER", "CSQ"),
         c("chrom", "pos", "ref", "alt", "mutect_filter", "vep_csq"))

fwrite(tumour_bam_read_counts, file = mutation_save_path)
