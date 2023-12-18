suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))

parser <- ArgumentParser()

parser$add_argument('--library_size_path',  nargs=1,
                    help='Path to library size',
                    required=TRUE)
parser$add_argument('--allele_read_count_path',  nargs=1,
                    help='Path to allele read count',
                    required=TRUE)
parser$add_argument('--sample_name',  nargs=1,
                    help='Sample name',
                    required=TRUE)
parser$add_argument('--reference_fasta_path' , nargs=1,
                    help='Path to the reference used in the alignment',
                    required=TRUE)
parser$add_argument('--out_file_prefix' , nargs=1,
                    help='Path to the reference used in the alignment',
                    required=TRUE)

args <- parser$parse_args()

library_size_path <- args$library_size_path
allele_read_count_path <- args$allele_read_count_path
sample_name <- args$sample_name
reference_fasta_path <- args$reference_fasta_path
out_file_prefix <- args$out_file_prefix

# get the read counts
allele_read_count <- fread(allele_read_count_path)

# get the library size
library_size <- as.numeric(readLines(library_size_path))
if(is.na(library_size)){
  stop("Library size is NA")
}
allele_read_count[, library_size := library_size]

# get the fasta reference and add the allele lengths
fasta <- read.fasta(reference_fasta_path)
for(line_idx in 1:nrow(allele_read_count)){
  
  a1 <- allele_read_count[line_idx]$allele1
  if(a1 != "" & !is.na(a1)){
    allele_read_count[line_idx, allele1_length := length(as.character(fasta[[a1]]))]
  }else{
    allele_read_count[line_idx, allele1_length := as.numeric(NA)]
  }
  
  a2 <- allele_read_count[line_idx]$allele2
  if(a2 != "" & !is.na(a2)){
    allele_read_count[line_idx, allele2_length := length(as.character(fasta[[a2]]))]
  }else{
    allele_read_count[line_idx, allele2_length := as.numeric(NA)]
  }
  
}

# add the fraction of reads that map uniquely to an allele
allele_read_count[,total_read_count := sum(c(allele1_n_reads_mapping_uniquely,
                                             allele2_n_reads_mapping_uniquely,
                                             n_multi_allele_reads), na.rm = TRUE), by = "gene"]
allele_read_count[,unique_read_count := sum(c(allele1_n_reads_mapping_uniquely,
                                              allele2_n_reads_mapping_uniquely), na.rm = TRUE), by = "gene"]

allele_read_count[total_read_count > 0, 
                  frac_mapping_uniquely := unique_read_count/total_read_count]
allele_read_count[total_read_count == 0, frac_mapping_uniquely := 0]

# add the fraction that map to multi gene
allele_read_count[,total_gene_read_count := sum(c(allele1_n_reads_mapping_uniquely,
                                                  allele2_n_reads_mapping_uniquely,
                                                  n_multi_allele_reads,
                                                  n_multi_gene_reads), na.rm = TRUE), by = "gene"]
allele_read_count[total_gene_read_count > 0, 
                  frac_mapping_multi_gene := n_multi_gene_reads/total_gene_read_count]
allele_read_count[total_gene_read_count == 0,
                  frac_mapping_multi_gene := 0]

# add in allele1 and allele2 frac
allele_read_count[homozygous == FALSE & unique_read_count > 0, allele1_frac := (allele1_n_reads_mapping_uniquely)/(allele1_n_reads_mapping_uniquely + allele2_n_reads_mapping_uniquely)]
allele_read_count[homozygous == FALSE & unique_read_count == 0, allele1_frac := 0]

allele_read_count[homozygous == FALSE & unique_read_count > 0, allele2_frac := (allele2_n_reads_mapping_uniquely)/(allele1_n_reads_mapping_uniquely + allele2_n_reads_mapping_uniquely)]
allele_read_count[homozygous == FALSE & unique_read_count == 0, allele2_frac := 0]

# add the updated allele1 count
allele_read_count[homozygous == FALSE, allele1_updated_count := allele1_n_reads_mapping_uniquely + allele1_frac*n_multi_allele_reads]
allele_read_count[homozygous == FALSE, allele2_updated_count := allele2_n_reads_mapping_uniquely + allele2_frac*n_multi_allele_reads]
allele_read_count[homozygous == TRUE & allele1 != "", allele1_updated_count := allele1_n_reads_mapping_uniquely]
allele_read_count[homozygous == TRUE & allele2 != "", allele2_updated_count := allele2_n_reads_mapping_uniquely]
                  
# add allele rpkm
allele_read_count[,allele1_length_kb := allele1_length/1000]
allele_read_count[,allele2_length_kb := allele2_length/1000]
allele_read_count[,library_size_scaled := library_size/1000000]

allele_read_count[homozygous == FALSE & allele1_updated_count > 0, allele1_rpkm := 
                    (allele1_updated_count/library_size_scaled)/allele1_length_kb]
allele_read_count[homozygous == FALSE & allele1_updated_count == 0, allele1_rpkm := 0]

allele_read_count[homozygous == FALSE & allele2_updated_count > 0, allele2_rpkm := 
                    (allele2_updated_count/library_size_scaled)/allele2_length_kb]
allele_read_count[homozygous == FALSE & allele2_updated_count == 0, allele2_rpkm := 0]

allele_read_count[homozygous == TRUE & allele1 != "" & allele1_updated_count > 0, 
                  allele1_rpkm := (allele1_updated_count/library_size_scaled)/allele1_length_kb]
allele_read_count[homozygous == TRUE & allele1 != "" & allele1_updated_count == 0, 
                  allele1_rpkm := 0]

allele_read_count[homozygous == TRUE & allele2 != "" & allele2_updated_count > 0, 
                  allele2_rpkm := (allele2_updated_count/library_size_scaled)/allele2_length_kb]
allele_read_count[homozygous == TRUE & allele2 != "" & allele2_updated_count == 0, 
                  allele2_rpkm := 0]

# add gene rpkm
allele_read_count[homozygous == FALSE, gene_rpkm := allele1_rpkm + allele2_rpkm]
allele_read_count[homozygous == TRUE & allele1 != "", gene_rpkm := allele1_rpkm]
allele_read_count[homozygous == TRUE & allele2 != "", gene_rpkm := allele2_rpkm]

allele_read_count <- allele_read_count[,c("patient", "sample_name", "gene", "allele1", "allele2", "homozygous",
                                          "allele1_n_reads_mapping_uniquely", "allele2_n_reads_mapping_uniquely", 
                                          "frac_mapping_multi_gene",
                                          "frac_mapping_uniquely", "allele1_updated_count", "allele2_updated_count",
                                          "allele1_length_kb", "allele2_length_kb", "library_size_scaled",
                                          "allele1_rpkm", "allele2_rpkm", "gene_rpkm")]

fwrite(allele_read_count, file = paste0(out_file_prefix, "_rpkm.csv"))
