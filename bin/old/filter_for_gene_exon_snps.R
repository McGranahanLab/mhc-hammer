suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser()

parser$add_argument('--gtf_path',  nargs=1,
                    help="Path to patient's gtf file",
                    required=TRUE)
parser$add_argument('--allele1_snp_path', nargs=1,
                    help="Path to .csv file containing snp positions.",
                    required=TRUE)
parser$add_argument('--allele2_snp_path', nargs=1,
                    help="Path to .csv file containing snp positions.",
                    required=TRUE)

args <- parser$parse_args()

gtf_path <- args$gtf_path
allele1_snp_path <- args$allele1_snp_path
allele2_snp_path <- args$allele2_snp_path

gtf <- fread(gtf_path, header = FALSE)
gtf <- gtf[,c("V1", "V3", "V4", "V5")]
setnames(gtf, c("allele", "feature_type", "start", "stop"))

# make allele 1 snps
allele1_snps <- fread(allele1_snp_path)
setnames(allele1_snps, c("allele", "start", "stop"))

allele1_name <- unique(allele1_snps$allele)
if(length(allele1_name) != 1){
  stop("Allele should be length")
}

for(line_idx in 1:nrow(allele1_snps)){
  snp_pos <- allele1_snps[line_idx]$stop
  feature <- gtf[allele == allele1_name & 
                   start <= snp_pos &
                   stop >= snp_pos]$feature_type
  if(length(feature) != 1){
    stop("Should be a single feature for a given snp")
  }
  allele1_snps[line_idx, feature_type := feature]
}

# make allele 2 snps
allele2_snps <- fread(allele2_snp_path)
setnames(allele2_snps, c("allele", "start", "stop"))

allele2_name <- unique(allele2_snps$allele)
if(length(allele2_name) != 1){
  stop("Allele should be length")
}

for(line_idx in 1:nrow(allele2_snps)){
  snp_pos <- allele2_snps[line_idx]$stop
  feature <- gtf[allele == allele2_name & 
                   start <= snp_pos &
                   stop >= snp_pos]$feature_type
  if(length(feature) != 1){
    stop("Should be a single feature for a given snp")
  }
  allele2_snps[line_idx, feature_type := feature]
}

# get which snps are in exons in both alleles
allele1_exon_snps <- allele1_snps[,feature_type == "exon"]
allele2_exon_snps <- allele2_snps[,feature_type == "exon"]

allele1_allele2_exon_snps <- allele1_exon_snps & allele2_exon_snps

allele1_exon_snps_dt <- allele1_snps[allele1_allele2_exon_snps]
allele2_exon_snps_dt <- allele2_snps[allele1_allele2_exon_snps]

allele1_exon_snps_dt[,feature_type := NULL]
allele2_exon_snps_dt[,feature_type := NULL]

fwrite(allele1_exon_snps_dt, file = gsub(".snp_pos.bed$", ".exon_snp_pos.bed", allele1_snp_path), sep = "\t", col.names = FALSE)
fwrite(allele2_exon_snps_dt, file = gsub(".snp_pos.bed$", ".exon_snp_pos.bed", allele2_snp_path), sep = "\t", col.names = FALSE)
