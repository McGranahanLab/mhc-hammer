suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--allele1_snp_bed',  nargs=1,
                    help='Path to allele 1 snp positions',
                    required=TRUE)
parser$add_argument('--allele2_snp_bed',  nargs=1,
                    help='Path to allele 2 snp positions',
                    required=TRUE)

parser$add_argument('--allele1',  nargs=1,
                    help='Allele 1',
                    required=TRUE)
parser$add_argument('--allele2',  nargs=1,
                    help='Allele2',
                    required=TRUE)

parser$add_argument('--allele1_filtered_pos_bed',  nargs=1,
                    help='Path to allele1 germline filtered coverage csv',
                    required=TRUE)
parser$add_argument('--allele2_filtered_pos_bed',  nargs=1,
                    help='Path to allele2 germline filtered coverage csv',
                    required=TRUE)

parser$add_argument('--sample_name',  nargs=1,
                    help='sample name',
                    required=TRUE)

args <- parser$parse_args()

allele1_snp_bed <- args$allele1_snp_bed
allele2_snp_bed <- args$allele2_snp_bed

allele1 <- args$allele1
allele2 <- args$allele2

allele1_filtered_pos_bed <- args$allele1_filtered_pos_bed
allele2_filtered_pos_bed <- args$allele2_filtered_pos_bed

sample_name <- args$sample_name

allele1_filtered_pos <- fread(allele1_filtered_pos_bed, sep = "\t")
allele2_filtered_pos <- fread(allele2_filtered_pos_bed, sep = "\t")

allele1_snp <- fread(allele1_snp_bed, sep = "\t")
allele2_snp <- fread(allele2_snp_bed, sep = "\t")

if(nrow(allele1_filtered_pos) > 0 & nrow(allele2_filtered_pos) > 0 &
   nrow(allele1_snp) > 0 & nrow(allele2_snp) > 0){
  
  setnames(allele1_filtered_pos, c("chr", "start", "end"))
  setnames(allele2_filtered_pos, c("chr", "start", "end"))  
  
  setnames(allele1_snp, c("chr", "start", "end"))
  setnames(allele2_snp, c("chr", "start", "end"))
  
  allele1_snp_positions_pass_cov <- allele1_snp[,end] %in% allele1_filtered_pos[,end]
  allele2_snp_positions_pass_cov <- allele2_snp[,end] %in% allele2_filtered_pos[,end]
  
  both_alleles_snp_positions_pass_cov <- allele1_snp_positions_pass_cov & allele2_snp_positions_pass_cov
  
  if(sum(both_alleles_snp_positions_pass_cov) > 0){
    snp_filtered_allele1 <- allele1_snp[,end][both_alleles_snp_positions_pass_cov]  
    snp_filtered_allele1_bed <- data.table(chr = allele1,
                                                 start = snp_filtered_allele1 - 1,
                                                 end = snp_filtered_allele1)
    
    snp_filtered_allele2 <- allele2_snp[,end][both_alleles_snp_positions_pass_cov]
    snp_filtered_allele2_bed <- data.table(chr = allele2,
                                                 start = snp_filtered_allele2 - 1,
                                                 end = snp_filtered_allele2)
    
    
    
    
  }else{
    snp_filtered_allele1_bed <- data.table(chr = character(),
                                                 start = numeric(),
                                                 end = numeric())
    snp_filtered_allele2_bed <- data.table(chr = character(),
                                                 start = numeric(),
                                                 end = numeric())
  }
  
}else{
  
  snp_filtered_allele1_bed <- data.table(chr = character(),
                                               start = numeric(),
                                               end = numeric())
  
  snp_filtered_allele2_bed <- data.table(chr = character(),
                                               start = numeric(),
                                               end = numeric())
}

fwrite(snp_filtered_allele1_bed,	
       file = paste0(sample_name, ".", allele1, ".filtered_snp_positions.bed"), 	
       sep="\t", col.names = FALSE)

fwrite(snp_filtered_allele2_bed,	
       file = paste0(sample_name, ".", allele2, ".filtered_snp_positions.bed"), 	
       sep="\t", col.names = FALSE)
