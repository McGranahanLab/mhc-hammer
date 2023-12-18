suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser()

parser$add_argument('--gtf_path',  nargs=1,
                    help="Path to patient's gtf file",
                    required=TRUE)
parser$add_argument('--snp_path', nargs=1,
                    help="Path to .csv file containing snp positions.",
                    required=TRUE)

args <- parser$parse_args()

gtf_path <- args$gtf_path
snp_path <- args$snp_path

gtf <- fread(gtf_path, header = FALSE)
gtf <- gtf[,c("V1", "V3", "V4", "V5")]
setnames(gtf, c("allele", "feature_type", "start", "stop"))

snps <- fread(snp_path)
setnames(snps, c("allele", "start", "stop"))

snp_allele <- unique(snps$allele)
if(length(snp_allele) != 1){
  stop("Allele should be length w")
}

for(line_idx in 1:nrow(snps)){
  snp_pos <- snps[line_idx]$stop
  feature <- gtf[allele == snp_allele & 
                   start <= snp_pos &
                   stop >= snp_pos]$feature_type
  if(length(feature) != 1){
    stop("Should be a single feature for a given snp")
  }
  snps[line_idx, feature_type := feature]
}

snps <- snps[feature_type == "exon"]

snps[,feature_type := NULL]

fwrite(snps, file = gsub(".snp_pos.bed$", ".exon_snp_pos.bed", snp_path), sep = "\t", col.names = FALSE)
