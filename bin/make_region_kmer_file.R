library(data.table)
library(seqinr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('--region_sj_table',  nargs=1,
                    help='Path to patient splice junction table',
                    required=TRUE)
parser$add_argument('--hla_genome_fasta_path',  nargs=1,
                    help='Path to HLA genome fasta path',
                    required=TRUE)
parser$add_argument('--gtf',  nargs=1,
                    help='Path to HLA gtf',
                    required=TRUE)
parser$add_argument('--sample_id',  nargs=1,
                    help='Patient name',
                    required=TRUE)
parser$add_argument('--original_kmer',  nargs=1,
                    help='Kmer file',
                    required=TRUE)

args <- parser$parse_args()

region_sj_table <- args$region_sj_table
gtf_path <- args$gtf
genome_fasta_path <- args$hla_genome_fasta_path
sample_id <- args$sample_id
kmer_file_path <- args$original_kmer

kmer_length <- 14
# check the kmer is even
if(!(kmer_length %% 2) == 0) {
  cat("kmer_length should be even, changing input kmer_length of ", kmer_length, " to ", kmer_length + 1, "\n")
  kmer_length <- kmer_length + 1
} 

# get the fasta
genome_fasta <- read.fasta(genome_fasta_path, forceDNAtolower = FALSE)

# get the gtf
gtf <- fread(gtf_path)
setnames(gtf, c("seqname",	"source",	"type",	"feature_start",	"feature_end",
                "score",	"strand",	"frame",	"other_values"))
alleles <- unique(gtf$seqname)

known_sjs <- gtf[type == "intron", c("seqname", "feature_start", "feature_end", "type")]
known_sjs <- known_sjs[order(seqname, feature_start)]

star_sj_tab <- fread(region_sj_table)
star_sj_tab <- star_sj_tab[,c("V1", "V2", "V3")]
setnames(star_sj_tab,
         c("chr", "start", "end"))

if(nrow(star_sj_tab) == 0){
  stop("Splice junction table is empty")
}

# add in the known sjs to the star_sj_tab if they are missing
star_sj_tab <- unique(rbindlist(list(star_sj_tab,
                                     known_sjs[,c("seqname", "feature_start", "feature_end")]), use.names = FALSE))
star_sj_tab <- star_sj_tab[order(chr, start)]
star_sj_tab[,gene := gsub("(hla_.).*", "\\1", chr)]

all_kmers <- data.table()
for(line_idx in 1:nrow(star_sj_tab)){
  cat(line_idx, "/", nrow(star_sj_tab), "\n")
  allele <- star_sj_tab[line_idx]$chr
  sj_start <- star_sj_tab[line_idx]$start
  sj_end <- star_sj_tab[line_idx]$end
  
  line_known_sjs <- known_sjs[seqname == allele]
  
  # check if this is a known or novel sjs
  if(nrow(line_known_sjs[feature_start == sj_start & feature_end == sj_end]) == 0 ){
    
    # this is a novel sj
    # remove any known sj that overlaps the novel sj
    
    # novel sj start lies within known sj
    line_known_sjs[sj_start > feature_start &
                     sj_start < feature_end, overlap := TRUE]
    
    # novel sj end lies within known sj
    line_known_sjs[sj_end > feature_start &
                     sj_end < feature_end, overlap := TRUE]
    
    # known sj start lies within novel sj
    line_known_sjs[feature_start > sj_start &
                     feature_start < sj_end, overlap := TRUE]
    
    line_known_sjs[feature_end > sj_start &
                     feature_end < sj_end, overlap := TRUE]
    
    new_sj_tab <- line_known_sjs[is.na(overlap), c("seqname", "feature_start", "feature_end")]
    new_sj_tab <- rbindlist(list(new_sj_tab,
                                 data.table(seqname = allele,
                                            feature_start = sj_start,
                                            feature_end = sj_end)))
    new_sj_tab <- new_sj_tab[order(feature_start)]
    
  }else if(nrow(line_known_sjs[feature_start == sj_start & feature_end == sj_end]) == 1 ){
    
    new_sj_tab <- line_known_sjs[, c("seqname", "feature_start", "feature_end")]
    # known sj
    
  }else{
    stop("known sjs should be unique")
  }
  
  allele_seq <- toupper(as.character(genome_fasta[[allele]]))
  genome_dt <- data.table(pos = 1:length(allele_seq),
                          base = allele_seq)
  # mark which bases are not transcribed
  for(i in 1:nrow(new_sj_tab)){
    genome_dt[pos %in% new_sj_tab[i]$feature_start : new_sj_tab[i]$feature_end, transcribed := FALSE]
  }
  genome_dt[is.na(transcribed), transcribed := TRUE]
  
  # from the novel sj_start count kmer/2 backwards
  half_kmer <- kmer_length/2  
  kmer_first_half <- vector(mode = "character", length = half_kmer)
  genome_idx <- sj_start - 1 # this counts along the genome
  kmer_first_half_idx <- half_kmer # this counts along kmer_first_half
  while(any(kmer_first_half == "")){
    dt <- genome_dt[pos == genome_idx]
    if(nrow(dt) > 0){
      # genome_idx is in sequence range
      if(dt$transcribed){
        kmer_first_half[kmer_first_half_idx] <- dt$base
        kmer_first_half_idx <- kmer_first_half_idx - 1
      }
      genome_idx <- genome_idx -1
    }else{
      # sj_start_idx is out of sequence range
      kmer_first_half <- kmer_first_half[kmer_first_half != ""]
    }
  }
  
  # from the novel sj_end count kmer/2 forwards
  kmer_second_half <- vector(mode = "character", length = half_kmer)
  genome_idx <- sj_end + 1 # this counts along the genome
  kmer_second_half_idx <- 1 # this counts along kmer_first_half
  while(any(kmer_second_half == "")){
    dt <- genome_dt[pos == genome_idx]
    if(nrow(dt) > 0){
      # genome_idx is in sequence range
      if(dt$transcribed){
        kmer_second_half[kmer_second_half_idx] <- dt$base
        kmer_second_half_idx <- kmer_second_half_idx + 1
      }
      genome_idx <- genome_idx + 1
    }else{
      # sj_start_idx is out of sequence range
      kmer_second_half <- kmer_second_half[kmer_second_half != ""]
    }
  }
  
  complete_kmer <- c(kmer_first_half, kmer_second_half)
  rev_complete_kmer <- rev(complete_kmer)
  A_idx <- which(rev_complete_kmer == "A")
  C_idx <- which(rev_complete_kmer == "C")
  G_idx <- which(rev_complete_kmer == "G")
  T_idx <- which(rev_complete_kmer == "T")
  
  rev_complete_kmer[A_idx] <- "T"
  rev_complete_kmer[C_idx] <- "G"
  rev_complete_kmer[G_idx] <- "C"
  rev_complete_kmer[T_idx] <- "A"
  
  complete_kmer <- paste0(complete_kmer, collapse = "")
  rev_complete_kmer <- paste0(rev_complete_kmer, collapse = "")
  
  all_kmers <- rbindlist(list(all_kmers,
                              data.table(kmer = c(complete_kmer, rev_complete_kmer))))
  
}
all_kmers <- unique(all_kmers)

original_kmer_file <- fread(kmer_file_path, header = FALSE)
setnames(original_kmer_file, "V1", "kmer")

both_kmer_files <- unique(rbindlist(list(original_kmer_file, all_kmers)))

fwrite(both_kmer_files,
       file = paste0(sample_id, "_kmers.txt"),
       col.names = FALSE)
