library(data.table)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('--tumour_novel_sjs_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--tumour_known_sjs_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--normal_novel_sjs_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--normal_known_sjs_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--sample_name',  nargs=1,
                    required=TRUE)
parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

args <- parser$parse_args()

tumour_novel_sjs_path <- args$tumour_novel_sjs_path
tumour_known_sjs_path <- args$tumour_known_sjs_path
normal_novel_sjs_path <- args$normal_novel_sjs_path
normal_known_sjs_path <- args$normal_known_sjs_path
sample_name <- args$sample_name
scripts_dir <- args$scripts_dir
source(paste0(scripts_dir, "/alt_splicing_functions.R"))

cat("tumour_novel_sjs_path = ", tumour_novel_sjs_path, "\n")
cat("tumour_known_sjs_path = ", tumour_known_sjs_path, "\n")
cat("normal_novel_sjs_path = ", normal_novel_sjs_path, "\n")
cat("normal_known_sjs_path = ", normal_known_sjs_path, "\n")

# tumour_novel_sjs_path <- "/camp/project/proj-tracerx-lung/tctProjects/putticc/mhc_pipeline_mutations/LOHHLA/mhc_hammer_results_tracerx/LTX0038/hla_alternative_splicing/LTX0038_SU_T1-R2--32a64edcc47a_novel_splice_junctions.csv"
# tumour_known_sjs_path <- "/camp/project/proj-tracerx-lung/tctProjects/putticc/mhc_pipeline_mutations/LOHHLA/mhc_hammer_results_tracerx/LTX0038/hla_alternative_splicing/LTX0038_SU_T1-R2--32a64edcc47a_known_splice_junctions.csv"
# 
# normal_novel_sjs_path <- "/camp/project/proj-tracerx-lung/tctProjects/putticc/mhc_pipeline_mutations/LOHHLA/mhc_hammer_results_tracerx/LTX0038/hla_alternative_splicing/LTX0038_SU_N01_novel_splice_junctions.csv"
# normal_known_sjs_path <- "/camp/project/proj-tracerx-lung/tctProjects/putticc/mhc_pipeline_mutations/LOHHLA/mhc_hammer_results_tracerx/LTX0038/hla_alternative_splicing/LTX0038_SU_N01_known_splice_junctions.csv"

# read in the splice junction tbales
tumour_novel_sjs <- fread(tumour_novel_sjs_path)
tumour_known_sjs <- fread(tumour_known_sjs_path)
normal_novel_sjs <- fread(normal_novel_sjs_path)
normal_known_sjs <- fread(normal_known_sjs_path)

tumour_sample_name <- unique(tumour_novel_sjs$sample_name)
normal_sample_name <- unique(normal_novel_sjs$sample_name)

tumour_novel_sjs <- tumour_novel_sjs[,c("allele", "gene", "start", "end", "n_unique_reads", "canonical_sj_read_count",
                                        "intron_n_reads", "novel_transcript_proportion", "total_read_count",
                                        "novel_sj_cat", "canonical_sj_names", "sj_type",
                                        "framshift", "premature_stop", "exon_intron_name")]

normal_novel_sjs <- normal_novel_sjs[,c("allele", "gene", "start", "end", "n_unique_reads", "canonical_sj_read_count",
                                        "intron_n_reads", "novel_transcript_proportion", "total_read_count",
                                        "novel_sj_cat", "canonical_sj_names", "sj_type",
                                        "framshift", "premature_stop", "exon_intron_name")]

setnames(tumour_novel_sjs,
         old = c("n_unique_reads", "canonical_sj_read_count",
                 "intron_n_reads", "novel_transcript_proportion", "total_read_count"),
         new = c("tumour_n_unique_reads", "tumour_canonical_sj_read_count", 
                 "tumour_intron_n_reads", "tumour_novel_transcript_proportion", "tumour_total_read_count"))

setnames(normal_novel_sjs,
         old = c("n_unique_reads", "canonical_sj_read_count",
                 "intron_n_reads", "novel_transcript_proportion", "total_read_count"),
         new = c("normal_n_unique_reads", "normal_canonical_sj_read_count", 
                 "normal_intron_n_reads", "normal_novel_transcript_proportion", "normal_total_read_count"))

tumour_normal_sjs <- merge(tumour_novel_sjs, normal_novel_sjs,
                           by = c("allele", "gene", "start", "end", "novel_sj_cat", "sj_type",
                                  "canonical_sj_names", "framshift", "premature_stop", "exon_intron_name"), 
                           all = TRUE)

# if not in tumour or not in normal then not detected, and n_unique reads is zero
if(nrow(tumour_normal_sjs[is.na(normal_n_unique_reads) & is.na(tumour_n_unique_reads)]) > 0){
  stop("Why are they both zero?")
}

tumour_normal_sjs[is.na(normal_n_unique_reads), normal_n_unique_reads := 0]
tumour_normal_sjs[is.na(tumour_n_unique_reads), tumour_n_unique_reads := 0]
tumour_normal_sjs[is.na(normal_novel_transcript_proportion), normal_novel_transcript_proportion := 0]
tumour_normal_sjs[is.na(tumour_novel_transcript_proportion), tumour_novel_transcript_proportion := 0]

# for splice junctions that were only called in the tumour, add the canonical read count in the normal
# (and the other way around)
for(line_idx in 1:nrow(tumour_normal_sjs)){
  
  if(tumour_normal_sjs[line_idx]$novel_sj_cat == "other"){
    next()
  }
  
  # found in tumour but not normal
  if(is.na(tumour_normal_sjs[line_idx]$normal_canonical_sj_read_count)){
    
    sj_allele <- tumour_normal_sjs[line_idx]$allele
    intron_name <- tumour_normal_sjs[line_idx]$canonical_sj_names
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "single_exon_skip"){
      
      intron_left <- gsub(";.*", "", intron_name)
      intron_left_reads <- normal_known_sjs[known_feature_name == intron_left & allele == sj_allele]$n_unique_reads
      
      intron_right <- gsub(".*;", "", intron_name)
      intron_right_reads <- normal_known_sjs[known_feature_name == intron_right & allele == sj_allele]$n_unique_reads
      
      if(length(intron_left_reads) > 1 | length(intron_right_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_left_reads) == 0){
        intron_left_reads <- as.numeric(NA)
      }
      if(length(intron_right_reads) == 0){
        intron_right_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, normal_canonical_sj_read_count := as.numeric(mean(c(intron_left_reads, intron_right_reads), na.rm = TRUE))]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "end_exon_skip"){
      
      intron_reads <- normal_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, normal_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "start_exon_skip"){
      
      # get the number of reads to the left and right
      intron_reads <- normal_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, normal_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "start_intron_retained"){
      
      intron_reads <- normal_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, normal_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "end_intron_retained"){
      
      intron_reads <- normal_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, normal_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
  }
  
  # found in normal but not tumour
  if(is.na(tumour_normal_sjs[line_idx]$tumour_canonical_sj_read_count)){
    
    sj_allele <- tumour_normal_sjs[line_idx]$allele
    intron_name <- tumour_normal_sjs[line_idx]$canonical_sj_names
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "single_exon_skip"){
      
      intron_left <- gsub(";.*", "", intron_name)
      intron_left_reads <- tumour_known_sjs[known_feature_name == intron_left & allele == sj_allele]$n_unique_reads
      
      intron_right <- gsub(".*;", "", intron_name)
      intron_right_reads <- tumour_known_sjs[known_feature_name == intron_right & allele == sj_allele]$n_unique_reads
      
      if(length(intron_left_reads) > 1 | length(intron_right_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_left_reads) == 0){
        intron_left_reads <- as.numeric(NA)
      }
      if(length(intron_right_reads) == 0){
        intron_right_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, tumour_canonical_sj_read_count := as.numeric(mean(c(intron_left_reads, intron_right_reads), na.rm = TRUE))]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "end_exon_skip"){
      
      intron_reads <- tumour_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, tumour_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "start_exon_skip"){
      
      # get the number of reads to the left and right
      intron_reads <- tumour_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, tumour_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "start_intron_retained"){
      
      intron_reads <- tumour_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, tumour_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
    if(tumour_normal_sjs[line_idx]$novel_sj_cat == "end_intron_retained"){
      
      intron_reads <- tumour_known_sjs[known_feature_name == intron_name & allele == sj_allele]$n_unique_reads
      
      if(length(intron_reads) > 1){
        stop("Introns should be unique to alleles")
      }
      
      if(length(intron_reads) == 0){
        intron_reads <- as.numeric(NA)
      }
      
      tumour_normal_sjs[line_idx, tumour_canonical_sj_read_count := as.numeric(intron_reads)]
    }
    
  }
  
}


# add in the novel_transcript_proportion change
tumour_normal_sjs[, novel_transcript_proportion_change := (tumour_novel_transcript_proportion - normal_novel_transcript_proportion)]
tumour_normal_sjs[tumour_novel_transcript_proportion > normal_novel_transcript_proportion, novel_transcript_proportion_changev2 := (tumour_novel_transcript_proportion - normal_novel_transcript_proportion)/tumour_novel_transcript_proportion]
tumour_normal_sjs[tumour_novel_transcript_proportion <= normal_novel_transcript_proportion, novel_transcript_proportion_changev2 := -1*(normal_novel_transcript_proportion - tumour_novel_transcript_proportion)/normal_novel_transcript_proportion]

for(line_idx in 1:nrow(tumour_normal_sjs)){
  
  if(tumour_normal_sjs[line_idx]$novel_sj_cat == "other"){
    next()
  }
  
  tumour_no_as <- ceiling(tumour_normal_sjs[line_idx]$tumour_canonical_sj_read_count)
  tumour_as <- ceiling(tumour_normal_sjs[line_idx]$tumour_n_unique_reads)
  normal_no_as <- ceiling(tumour_normal_sjs[line_idx]$normal_canonical_sj_read_count)
  normal_as <- ceiling(tumour_normal_sjs[line_idx]$normal_n_unique_reads)
  
  if(any(is.na(c(tumour_no_as, tumour_as, normal_no_as, normal_as)))){
    next
  }
  
  m <- matrix(c(tumour_no_as, tumour_as, 
                normal_no_as, normal_as),
              nrow = 2, byrow = TRUE)
  ft <- fisher.test(m)
  
  tumour_normal_sjs[line_idx, fisher_pvalue := as.numeric(ft$p.value)]
  tumour_normal_sjs[line_idx, fisher_odds_ratio := as.numeric(ft$estimate)]
  tumour_normal_sjs[line_idx, fisher_ci_lower := as.numeric(ft$conf.int[1])]
  tumour_normal_sjs[line_idx, fisher_ci_upper := as.numeric(ft$conf.int[2])]
  
}

tumour_normal_sjs[,tumour_sample_name := tumour_sample_name]
tumour_normal_sjs[,normal_sample_name := normal_sample_name]

fwrite(tumour_normal_sjs, file = paste0(sample_name, "_tumour_normal_splice_junctions.csv"))

