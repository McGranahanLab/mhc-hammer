suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--snp_reads_overlap_bed',  nargs=1,
                    help='Bed file with snp/reads overlap',
                    required=TRUE)
parser$add_argument('--out_csv',  nargs=1,
                    help='Name for output coverage file where reads can only count towards a snp once',
                    required=TRUE)

args <- parser$parse_args()

snp_reads_overlap_bed <- args$snp_reads_overlap_bed
out_csv <- args$out_csv

snp_reads_overlap <- fread(snp_reads_overlap_bed, sep = "\t")
if(nrow(snp_reads_overlap) > 0){
  setnames(snp_reads_overlap, c("bed_chr", "bed_start", "bed_end",
                                "bam_chr", "bam_start", "bam_end", "read"))
  
  # get all positions so we can add back zero depth positions later
  # (positions can have zero depth either because no reads overlap them,
  # or because the only reads that overlap them are counted for other snps)
  all_snp_positions <- unique(snp_reads_overlap$bed_end)
  
  # remove null features
  snp_reads_overlap <- snp_reads_overlap[bam_chr != "."]
  
  if(nrow(snp_reads_overlap) > 0){
    
    # remove lines that have duplicated reads
    snp_reads_overlap[, dup_id := 1:.N, by = "read"]
    snp_reads_overlap_no_dupes <- snp_reads_overlap[dup_id == 1]
    
    # now get the depth at the SNPs
    coverage <- snp_reads_overlap_no_dupes[,.N,by = "bed_end"]
    setnames(coverage, c("pos", "depth"))
    
    # snp positions that dont have any reads (i.e. dup_id > 1 for all reads) need to be 
    # added back with depth zero
    missing_positions <- unique(all_snp_positions[!all_snp_positions %in% coverage$pos])
    
    if(length(missing_positions) > 0){
      tmp_dt <- data.table(pos = missing_positions,
                           depth = 0)
      coverage <- rbindlist(list(coverage, tmp_dt))
      coverage <- coverage[order(pos)]
    }
    
    coverage[, depth := depth + 1]
  }else{
    coverage <- data.table(pos = as.numeric(),
                           depth = as.numeric())
  }
}else{
  coverage <- data.table(pos = as.numeric(),
                         depth = as.numeric())
}

fwrite(coverage, file = out_csv)



