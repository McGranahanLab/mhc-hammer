suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--mpileup_path',  nargs=1,
                    help='Path to mpileup file',
                    required=TRUE)
parser$add_argument('--hla_ref_path',  nargs=1,
                    help='Path to HLA reference',
                    required=TRUE)
parser$add_argument('--allele',  nargs=1,
                    help='HLA allele',
                    required=TRUE)
parser$add_argument('--out_csv_name',  nargs=1,
                    help='Path to output csv',
                    required=TRUE)

args <- parser$parse_args()

mpileup_path <- args$mpileup_path
hla_ref_path <- args$hla_ref_path
allele <- args$allele
out_csv_name <- args$out_csv_name

get_mpileup <- function(mpileup_path){
  
  if(!file.exists(mpileup_path)){
    stop("Trying to read file that does not exist: ", mpileup_path)
  }
  
  allele_mpileup <- fread(mpileup_path, sep = "\t")
  if(nrow(allele_mpileup) == 0 & ncol(allele_mpileup) == 0){
    cat("mpileup empty\n")
  }else{
    setnames(allele_mpileup, c("allele", "pos", "base", "depth", "read_bases", "base_qualities"))
    allele_mpileup <- allele_mpileup[,c("allele", "pos", "base", "depth")]
  }
  
  return(allele_mpileup)
  
}

hla_ref   <- read.fasta(hla_ref_path)
allele_length <- length(hla_ref[[allele]])

mpileup <- get_mpileup(mpileup_path = mpileup_path)
if(ncol(mpileup) == 0 & nrow(mpileup) == 0){
  coverage_dt <- data.table(pos = 1:allele_length,
                            depth = 0)
}else{
  # double check we have all positions
  # some versions of samtools dont print zero depth positions
  missing_positions <- (1:allele_length)[!(1:allele_length %in% mpileup$pos)]
  
  if(length(missing_positions) > 0){
    missing_positions_dt <- data.table(pos = missing_positions, depth = 0)
    
    coverage_dt <- rbindlist(list(mpileup[,c("pos", "depth")],
                                  missing_positions_dt), use.names = TRUE)
  }else{
    coverage_dt <- mpileup[,c("pos", "depth")]
  }
  
  coverage_dt <- coverage_dt[order(pos)]
}
fwrite(coverage_dt, out_csv_name)




