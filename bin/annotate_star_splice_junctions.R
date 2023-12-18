suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(argparse))

########## code ###########
parser <- ArgumentParser()

parser$add_argument('--sj_tab_path',  nargs=1,
                    help='Path to patient splice junction table',
                    required=TRUE)
parser$add_argument('--hla_genome_fasta_path',  nargs=1,
                    help='Path to HLA genome fasta path',
                    required=TRUE)
parser$add_argument('--gtf_path',  nargs=1,
                    help='Path to HLA gtf',
                    required=TRUE)
parser$add_argument('--sample_id',  nargs=1,
                    help='Sample ID',
                    required=TRUE)
parser$add_argument('--codon_table_path',  nargs=1,
                    help='Sample ID',
                    required=TRUE)
parser$add_argument('--scripts_dir',  nargs=1,
                    help='Path to project scripts dir',
                    required=TRUE)

args <- parser$parse_args()

sj_tab_path <- args$sj_tab_path
gtf_path <- args$gtf_path
hla_genome_fasta_path <- args$hla_genome_fasta_path
sample_id <- args$sample_id
codon_table_path <- args$codon_table_path

scripts_dir <- args$scripts_dir
source(paste0(scripts_dir, "/alt_splicing_functions.R"))

cat("sj_tab_path = ", sj_tab_path, "\n")
cat("gtf_path = ", gtf_path, "\n")
cat("hla_genome_fasta_path = ", hla_genome_fasta_path, "\n")
cat("sample_id = ", sample_id, "\n")
cat("codon_table_path = ", codon_table_path, "\n")

star_sj_tab <- fread(sj_tab_path)
setnames(star_sj_tab, 
         c("allele", "start", "end", "strand", "intron_motif", "annotation",
           "n_unique_reads", "n_multtimap_reads", "max_overhang"))
star_sj_tab[, gene := gsub("(hla_.).*", "\\1", allele)]

# get the GTF files 
gtf <- fread(gtf_path)
setnames(gtf, c("seqname",	"source",	"type",	"feature_start",	"feature_end",	
                "score",	"strand",	"frame",	"other_values"))
gtf[, feature_name := gsub('.*feature_name "(.*?)".*', "\\1", other_values)]

# create a table of known splice sites
known_splice_sites <- gtf[type == "intron", c("seqname", "feature_start", "feature_end", "feature_name")]
setnames(known_splice_sites, old = "feature_name", new = "known_feature_name")
known_splice_sites[, splice_site_type := "known"]

# add known splice junctions to splice junction table
star_sj_tab <- merge(star_sj_tab, known_splice_sites, 
                     by.x = c("allele", "start", "end"),
                     by.y = c("seqname", "feature_start", "feature_end"),
                     all = TRUE)

# splice junctions not in gtf are novel
star_sj_tab[is.na(splice_site_type), splice_site_type :="novel"]

# splice junctions with NA unique reads have zero unique reads
star_sj_tab[is.na(n_unique_reads), n_unique_reads := 0]
star_sj_tab[is.na(n_multtimap_reads), n_multtimap_reads := 0]

# get novel and known table
known_sjs <- star_sj_tab[splice_site_type == "known"]
novel_sjs <- star_sj_tab[splice_site_type == "novel"]
novel_sjs[,known_feature_name := NULL]

# now annotate each novel site
for(line_idx in 1:nrow(novel_sjs)){
  
  sj_start <- novel_sjs[line_idx]$start
  sj_end <- novel_sjs[line_idx]$end
  sj_chr <- novel_sjs[line_idx]$allele
  gene <- novel_sjs[line_idx]$gene
  
  if(!gene %in% c("hla_a", "hla_b", "hla_c")){
    stop("Gene should be hla_a, hla_b or hla_c")
  }
  
  # check if splice junction lies completely within a single feature
  sj_within_feature <- check_sj_within_feature(sj_chr, sj_start, sj_end, gtf)
  novel_sjs[line_idx, sj_in_feature := sj_within_feature$sj_in_feature]
  novel_sjs[line_idx, sj_in_feature_name := sj_within_feature$sj_in_feature_name]
  
  # check for full exon skiping
  exon_skip <- check_exon_skip(sj_chr, sj_start, sj_end, gtf)
  novel_sjs[line_idx, full_exon_skip := exon_skip$full_exon_skip]
  novel_sjs[line_idx, n_exon_skipped := exon_skip$n_exon_skipped]
  novel_sjs[line_idx, exon_skipped_names := exon_skip$exon_skipped_names]   
  
  # check if the end of an exon is skipped
  if(gene == "hla_a"){
    exon_end_skip <- check_exon_end_skip_forward_strand(sj_chr, sj_start, gtf)
  }else if(gene %in% c("hla_b",  "hla_c")){
    exon_end_skip <- check_exon_end_skip_reverse_strand(sj_chr, sj_end, gtf)
  }
  novel_sjs[line_idx, end_exon_skip := exon_end_skip$end_exon_skip]
  novel_sjs[line_idx, end_exon_skipped_name := exon_end_skip$end_exon_skipped_name]
  novel_sjs[line_idx, end_exon_skipped_start := exon_end_skip$end_exon_skipped_start]
  novel_sjs[line_idx, end_exon_skipped_end := exon_end_skip$end_exon_skipped_end]
  
  # check if exon start skipped
  if(gene == "hla_a"){
    exon_start_skip <- check_exon_start_skip_forward_strand(sj_chr, sj_end, gtf)
  }else if(gene %in% c("hla_b",  "hla_c")){
    exon_start_skip <- check_exon_start_skip_reverse_strand(sj_chr, sj_start, gtf)
  }
  novel_sjs[line_idx, start_exon_skip := exon_start_skip$start_exon_skip]
  novel_sjs[line_idx, start_exon_skipped_name := exon_start_skip$start_exon_skipped_name]
  novel_sjs[line_idx, start_exon_skipped_start := exon_start_skip$start_exon_skipped_start]
  novel_sjs[line_idx, start_exon_skipped_end := exon_start_skip$start_exon_skipped_end]
  
  # check for partial intron start retention
  if(gene == "hla_a"){
    intron_start_retained <- check_intron_start_retained_forward_strand(sj_chr, sj_start, gtf)
  }else if(gene %in% c("hla_b",  "hla_c")){
    intron_start_retained <- check_intron_start_retained_reverse_strand(sj_chr, sj_end, gtf)
  }
  novel_sjs[line_idx, start_intron_retained := intron_start_retained$start_intron_retained]
  novel_sjs[line_idx, start_intron_retained_name := intron_start_retained$start_intron_retained_name]
  novel_sjs[line_idx, start_intron_retained_start := intron_start_retained$start_intron_retained_start]
  novel_sjs[line_idx, start_intron_retained_end := intron_start_retained$start_intron_retained_end]
  
  # check for partial intron end retention
  if(gene == "hla_a"){
    intron_end_retained <- check_intron_end_retained_forward_strand(sj_chr, sj_end, gtf)
  }else if(gene %in% c("hla_b",  "hla_c")){
    intron_end_retained <- check_intron_end_retained_reverse_strand(sj_chr, sj_start, gtf)
  }
  novel_sjs[line_idx, end_intron_retained := intron_end_retained$end_intron_retained]
  novel_sjs[line_idx, end_intron_retained_name := intron_end_retained$end_intron_retained_name]
  novel_sjs[line_idx, end_intron_retained_start := intron_end_retained$end_intron_retained_start]
  novel_sjs[line_idx, end_intron_retained_end := intron_end_retained$end_intron_retained_end]
}

# add in overall categories
novel_sjs[ sj_in_feature == FALSE &
             full_exon_skip == TRUE &
             end_exon_skip == FALSE &
             start_exon_skip == FALSE &
             start_intron_retained == FALSE &
             end_intron_retained == FALSE &
             n_exon_skipped == 1, 
           novel_sj_cat := "single_exon_skip"]

novel_sjs[sj_in_feature == FALSE &
            full_exon_skip == FALSE &
            end_exon_skip == TRUE &
            start_exon_skip == FALSE &
            start_intron_retained == FALSE &
            end_intron_retained == FALSE, 
          novel_sj_cat := "end_exon_skip"]

novel_sjs[sj_in_feature == FALSE &
            full_exon_skip == FALSE &
            end_exon_skip == FALSE &
            start_exon_skip == TRUE &
            start_intron_retained == FALSE &
            end_intron_retained == FALSE, 
          novel_sj_cat := "start_exon_skip"]

novel_sjs[sj_in_feature == FALSE &
            full_exon_skip == FALSE &
            end_exon_skip == FALSE &
            start_exon_skip == FALSE &
            start_intron_retained == TRUE &
            end_intron_retained == FALSE, 
          novel_sj_cat := "start_intron_retained"]

novel_sjs[sj_in_feature == FALSE &
            full_exon_skip == FALSE &
            end_exon_skip == FALSE &
            start_exon_skip == FALSE &
            start_intron_retained == FALSE &
            end_intron_retained == TRUE, 
          novel_sj_cat := "end_intron_retained"]

novel_sjs[is.na(novel_sj_cat), 
          novel_sj_cat := "other"]

# add the consequence
hla_genome_fasta <- read.fasta(hla_genome_fasta_path, forceDNAtolower = FALSE)
codon_aa_tab <- fread(codon_table_path, header = TRUE)

for(line_idx in 1:nrow(novel_sjs)){
  
  if(novel_sjs[line_idx]$novel_sj_cat == "other"){
    next()
  }
  
  allele <- novel_sjs[line_idx]$allele
  gene <- novel_sjs[line_idx]$gene
  
  novel_start <- novel_sjs[line_idx]$start
  novel_end <- novel_sjs[line_idx]$end
  
  allele_gtf <- gtf[seqname == allele]
  genome_seq <- toupper(as.character(hla_genome_fasta[[allele]]))
  
  # check the coordinates match
  if(max(allele_gtf$feature_end) != length(genome_seq)){
    stop("The coordinates of the GTF and fasta dont match")
  }
  
  # add the  sequence to the gtf
  for(gtf_line in 1:nrow(allele_gtf)){
    start <- allele_gtf[gtf_line]$feature_start
    end <- allele_gtf[gtf_line]$feature_end
    allele_gtf[gtf_line, seq := paste0(genome_seq[start:end], collapse = "")]
  }
  
  # convert original sequence to protein
  original_seq <- toupper(unlist(strsplit(allele_gtf$seq, "")))
  original_exon_seq <- toupper(unlist(strsplit(allele_gtf[type == "exon"]$seq, "")))
  
  if(gene %in% c("hla_b", "hla_c")){
    original_exon_seq_dt <- data.table(pos = 1:length(original_exon_seq),
                                       base = original_exon_seq)
    original_exon_seq_dt[base == "G", new_base := "C"]
    original_exon_seq_dt[base == "C", new_base := "G"]
    original_exon_seq_dt[base == "A", new_base := "T"]
    original_exon_seq_dt[base == "T", new_base := "A"]
    original_exon_seq <- rev(original_exon_seq_dt$new_base)
  }
  
  original_codon_dt <- get_codon_tab(original_exon_seq, codon_aa_tab)
  original_protein <- original_codon_dt$amino_acid
  
  if(original_protein[length(original_protein)] != "*"){
    warning("Protein doesnt end in stop codon")
  }
  
  # get the new protein that results from the alternate splicing
  # if the novel sj involves intron retention we need to update the coordinates
  # otherwise we just add the coordinates
  novel_allele_gtf <- copy(allele_gtf)
  
  novel_allele_gtf[, new_feature_start := as.numeric(NA)]
  novel_allele_gtf[, new_feature_end := as.numeric(NA)]
  
  if(novel_sjs[line_idx]$novel_sj_cat %in% c("end_intron_retained", "start_intron_retained") ){
    
    if(novel_sjs[line_idx]$novel_sj_cat == "start_intron_retained" & gene == "hla_a"){
      # update the intron start to be the splice start
      novel_allele_gtf[type == "intron" & 
                         ( novel_start >= feature_start &
                             novel_start <= feature_end), 
                       new_feature_start :=  novel_start]
      
    }
    
    if(novel_sjs[line_idx]$novel_sj_cat == "end_intron_retained" & gene == "hla_a"){
      
      novel_allele_gtf[type == "intron" & 
                         ( novel_end >= feature_start &
                             novel_end <= feature_end), 
                       new_feature_end :=  novel_end]
      
    }
    
    if(novel_sjs[line_idx]$novel_sj_cat == "start_intron_retained" & gene %in% c("hla_b", "hla_c")){
      
      novel_allele_gtf[type == "intron" & 
                         ( novel_end >= feature_start &
                             novel_end <= feature_end), 
                       new_feature_end :=  novel_end]
      
    }
    
    if(novel_sjs[line_idx]$novel_sj_cat == "end_intron_retained" & gene %in% c("hla_b", "hla_c")){
      
      novel_allele_gtf[type == "intron" & 
                         ( novel_start >= feature_start &
                             novel_start <= feature_end), 
                       new_feature_start :=  novel_start]
      
    }
    
    # check we have updated the coordinates of a single intron
    if( ( nrow(novel_allele_gtf[!is.na(new_feature_start)]) +
          nrow(novel_allele_gtf[!is.na(new_feature_end)]) ) != 1){
      stop("Looks like the coordinates of more than one intron has been updated") 
    }
    
  }else{
    # if intron not retained just add new sj as intron
    # this will remove the bases in novel splice junction
    novel_allele_gtf <- rbindlist(list(novel_allele_gtf, 
                                       data.table(seqname = allele,
                                                  type = "intron",
                                                  feature_start = novel_start,
                                                  feature_end = novel_end)), 
                                  use.names = TRUE, fill = TRUE)
  }
  
  novel_allele_gtf[is.na(new_feature_start), new_feature_start := feature_start]
  novel_allele_gtf[is.na(new_feature_end), new_feature_end := feature_end]
  
  # make a dt of all the bases
  full_seq_dt <- data.table(pos = 1:length(original_seq),
                            base = original_seq)
  # if it falls in a UTR or intron it is not translated
  for(idx in 1:nrow(novel_allele_gtf)){
    if(novel_allele_gtf[idx]$type %in% c("intron", "5UTR", "3UTR")){
      full_seq_dt[pos >= novel_allele_gtf[idx]$new_feature_start &
                    pos <= novel_allele_gtf[idx]$new_feature_end , translated := FALSE]
    }
  }  
  
  full_seq_dt[is.na(translated), translated := TRUE]
  
  # quick check
  translated_count <- full_seq_dt[,.N,by = "translated"]
  if(nrow(translated_count) != 2){
    stop("Either all translated or none translated?")
  }
  
  new_exon_seq <- full_seq_dt[translated == TRUE]$base
  if(gene %in% c("hla_b", "hla_c")){
    new_exon_seq_dt <- data.table(pos = 1:length(new_exon_seq),
                                  base = new_exon_seq)
    new_exon_seq_dt[base == "G", new_base := "C"]
    new_exon_seq_dt[base == "C", new_base := "G"]
    new_exon_seq_dt[base == "A", new_base := "T"]
    new_exon_seq_dt[base == "T", new_base := "A"]
    new_exon_seq <- rev(new_exon_seq_dt$new_base)
  }
  
  new_codon_dt <- get_codon_tab(toupper(new_exon_seq), codon_aa_tab)
  new_protein <- new_codon_dt$amino_acid
  
  # check the difference between them
  seq_dist <- Biostrings::stringDist(c(paste0(new_protein, collapse = ""), 
                                       paste0(original_protein, collapse = "")))
  lev_dist <- as.numeric(seq_dist)
  
  # is a novel stop codon introduced
  ptc <- FALSE
  if("*" %in% new_protein[1:(length(new_protein) - 1)]){
    ptc <- TRUE  
    # get the position of the ptc
  }
  novel_sjs[line_idx, premature_stop := ptc]
  novel_sjs[line_idx, original_seq := paste0(original_exon_seq, collapse = "")]
  novel_sjs[line_idx, new_seq := paste0(toupper(full_seq_dt[translated == TRUE]$base), collapse = "")]
  novel_sjs[line_idx, original_amino_acid := paste0(original_protein, collapse = "")]
  novel_sjs[line_idx, new_amino_acid := paste0(new_protein, collapse = "")]
  novel_sjs[line_idx, protein_lev_dist := lev_dist]
  
}

# add the number of reads that support the canonical splice junction
for(line_idx in 1:nrow(novel_sjs)){
  
  if(novel_sjs[line_idx]$novel_sj_cat == "other"){
    next()
  }
  
  sj_allele <- novel_sjs[line_idx]$allele
  
  if(novel_sjs[line_idx]$novel_sj_cat == "single_exon_skip"){
    
    # get the number of reads to the left and right
    skipped_exon <- as.numeric(gsub("exon_", "", novel_sjs[line_idx]$exon_skipped_names))
    
    intron_left <- paste0("intron_", skipped_exon - 1)
    intron_left_reads <- known_sjs[known_feature_name == intron_left & allele == sj_allele]$n_unique_reads
    
    intron_right <- paste0("intron_", skipped_exon)
    intron_right_reads <- known_sjs[known_feature_name == intron_right & allele == sj_allele]$n_unique_reads
    
    if(length(intron_left_reads) > 1 | length(intron_right_reads) > 1){
      stop("Introns should be unique to alleles")
    }
    
    if(length(intron_left_reads) == 0){
      intron_left_reads <- as.numeric(NA)
    }
    if(length(intron_right_reads) == 0){
      intron_right_reads <- as.numeric(NA)
    }
    
    novel_sjs[line_idx, canonical_sj_names := paste0(intron_left, ";", intron_right)]  
    novel_sjs[line_idx, canonical_sj_read_count := as.numeric(mean(c(intron_left_reads, intron_right_reads), na.rm = TRUE))]
  }
  
  if(novel_sjs[line_idx]$novel_sj_cat == "end_exon_skip"){
    
    # get the number of reads to the left and right
    skipped_exon <- as.numeric(gsub("exon_", "", novel_sjs[line_idx]$end_exon_skipped_name))
    intron <- paste0("intron_", skipped_exon)
    intron_reads <- known_sjs[known_feature_name == intron & allele == sj_allele]$n_unique_reads
    
    if(length(intron_reads) > 1){
      stop("Introns should be unique to alleles")
    }
    
    if(length(intron_reads) == 0){
      intron_reads <- as.numeric(NA)
    }
    
    novel_sjs[line_idx, canonical_sj_names := intron]  
    novel_sjs[line_idx, intron_n_reads := intron_reads]
    novel_sjs[line_idx, canonical_sj_read_count := as.numeric(intron_reads)]
  }
  
  if(novel_sjs[line_idx]$novel_sj_cat == "start_exon_skip"){
    
    # get the number of reads to the left and right
    skipped_exon <- as.numeric(gsub("exon_", "", novel_sjs[line_idx]$start_exon_skipped_name))
    intron <- paste0("intron_", (skipped_exon - 1))
    intron_reads <- known_sjs[known_feature_name == intron & allele == sj_allele]$n_unique_reads
    
    if(length(intron_reads) > 1){
      stop("Introns should be unique to alleles")
    }
    
    if(length(intron_reads) == 0){
      intron_reads <- as.numeric(NA)
    }
    
    novel_sjs[line_idx, canonical_sj_names := intron]  
    novel_sjs[line_idx, intron_n_reads := intron_reads]
    novel_sjs[line_idx, canonical_sj_read_count := as.numeric(intron_reads)]
  }
  
  if(novel_sjs[line_idx]$novel_sj_cat == "start_intron_retained"){
    
    # get the number of reads to the left and right
    retained_intron <- as.numeric(gsub("intron_", "", novel_sjs[line_idx]$start_intron_retained_name))
    intron <- paste0("intron_", retained_intron)
    intron_reads <- known_sjs[known_feature_name == intron & allele == sj_allele]$n_unique_reads
    
    if(length(intron_reads) > 1){
      stop("Introns should be unique to alleles")
    }
    
    if(length(intron_reads) == 0){
      intron_reads <- as.numeric(NA)
    }
    
    novel_sjs[line_idx, canonical_sj_names := intron]  
    novel_sjs[line_idx, intron_n_reads := intron_reads]
    novel_sjs[line_idx, canonical_sj_read_count := as.numeric(intron_reads)]
  }
  
  if(novel_sjs[line_idx]$novel_sj_cat == "end_intron_retained"){
    
    # get the number of reads to the left and right
    retained_intron <- as.numeric(gsub("intron_", "", novel_sjs[line_idx]$end_intron_retained_name))
    intron <- paste0("intron_", retained_intron)
    intron_reads <- known_sjs[known_feature_name == intron & allele == sj_allele]$n_unique_reads
    
    if(length(intron_reads) > 1){
      stop("Introns should be unique to alleles")
    }
    
    if(length(intron_reads) == 0){
      intron_reads <- as.numeric(NA)
    }
    
    novel_sjs[line_idx, canonical_sj_names := intron]  
    novel_sjs[line_idx, intron_n_reads := intron_reads]
    novel_sjs[line_idx, canonical_sj_read_count := as.numeric(intron_reads)]
  }
}

# add the novel to canonical ratio
novel_sjs[, total_read_count := canonical_sj_read_count + n_unique_reads]
novel_sjs[, ratio := n_unique_reads/total_read_count]

# add a broad category 
novel_sjs[novel_sj_cat %in% c("end_intron_retained", "start_intron_retained"),
          sj_type := "partial_intron_retention"]
novel_sjs[novel_sj_cat %in% c("start_exon_skip", "end_exon_skip"),
          sj_type := "partial_exon_skip"]
novel_sjs[novel_sj_cat == "single_exon_skip",
          sj_type := "complete_exon_skip"]
novel_sjs[novel_sj_cat == "other",
          sj_type := "other"]

# add in a single column with the exon/intron involved
novel_sjs[novel_sj_cat == "start_exon_skip", exon_intron_name := start_exon_skipped_name]
novel_sjs[novel_sj_cat == "end_exon_skip", exon_intron_name := end_exon_skipped_name]

novel_sjs[novel_sj_cat == "start_intron_retained", exon_intron_name := start_intron_retained_name]
novel_sjs[novel_sj_cat == "end_intron_retained", exon_intron_name := end_intron_retained_name]

novel_sjs[novel_sj_cat == "single_exon_skip", exon_intron_name := exon_skipped_names]

# add in if there is a frameshift event
novel_sjs[new_seq != "", new_seq_length := nchar(new_seq)]
novel_sjs[!is.na(new_seq_length) & new_seq_length %% 3 == 0, framshift := FALSE]
novel_sjs[!is.na(new_seq_length) & new_seq_length %% 3 != 0, framshift := TRUE]
novel_sjs[,new_seq_length := NULL]

fwrite(novel_sjs, file = paste0(sample_id, "_novel_splice_junctions.csv"))
fwrite(known_sjs, file = paste0(sample_id, "_known_splice_junctions.csv"))


