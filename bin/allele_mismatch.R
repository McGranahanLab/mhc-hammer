#!/usr/bin/env 

##### Inputs ####

suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

##### Functions #####

getMisMatchPositionsPairwiseAlignment <- function(alignment){
  
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  
  sequence1 <- unlist(strsplit(as.character(seq1aln),split=""))
  sequence2 <- unlist(strsplit(as.character(seq2aln),split=""))
  
  k <- 1
  seq1Positions   <- c()
  for (char in sequence1){
    if(char%in%c('C','G','A','T')){
      seq1Positions <- c(seq1Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-')){
      # stop()
      seq1Positions <- c(seq1Positions,k)
      next;
      #
    }
  }
  
  
  k <- 1
  seq2Positions   <- c()
  for (char in sequence2){
    if(char%in%c('C','G','A','T')){
      seq2Positions <- c(seq2Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-')){
      seq2Positions <- c(seq2Positions,k)
      next;
      #
    }
  }
  
  diffSeq1 <-  seq1Positions[sequence1!=sequence2]
  diffSeq2 <-  seq2Positions[sequence1!=sequence2]
  
  diffType1 <- rep(1,length(diffSeq1))
  diffType1[which(sequence1[sequence1!=sequence2]%in%'-')] <- 2
  
  diffType2 <- rep(1,length(diffSeq2))
  diffType2[which(sequence2[sequence2!=sequence1]%in%'-')] <- 2
  
  
  out <- list()
  out$diffSeq1 <- diffSeq1
  out$diffSeq2 <- diffSeq2
  out$diffType1 <- diffType1
  out$diffType2 <- diffType2
  return(out)
  
}

get_hla_in_fasta <- function(allele, hla_fasta_alleles){
  
  if(!allele %in% hla_fasta_alleles){
    alt_hla_allele <- grep(pattern = allele, x = hla_fasta_alleles, value = TRUE)[1]
    return(alt_hla_allele)
  }else{
    return(allele)
  }
}

parser <- ArgumentParser()

parser$add_argument('--patient_id',  nargs=1,
                    help='Patient ID.',
                    required=TRUE)
parser$add_argument('--patient_hla_alleles', nargs=1,
                    help="Path to .csv file containing patient's HLA calls.",
                    required=TRUE)
parser$add_argument('--mhc_fasta', nargs=1,
                    help="Path to complete HLA allele fasta file",
                    required=TRUE)
parser$add_argument('--genome_or_transcriptome', nargs=1,
                    help="Are the mismatch positions based on the genome or transcriptome",
                    required=TRUE)

args <- parser$parse_args()

patient_id <- args$patient_id
patient_hla_alleles_path <- args$patient_hla_alleles
mhc_fasta_path <- args$mhc_fasta
genome_transcriptome <- args$genome_or_transcriptome

if(!file.exists(patient_hla_alleles_path)){
  stop("File does not exist: ", patient_hla_alleles_path)
}

if(!file.exists(mhc_fasta_path)){
  stop("File does not exist: ", mhc_fasta_path)
}

##### Get patients HLA calls #####

cat("Reading HLA calls input: ", patient_hla_alleles_path, "\n")

patient_alleles_dt <- fread(patient_hla_alleles_path, header = FALSE)
setnames(patient_alleles_dt, c("gene", "allele1", "allele2"))
patient_alleles_dt <- patient_alleles_dt[allele1 != "not typed" & allele2 != "not typed"]

# add the patient to the table
patient_alleles_dt[, patient := patient_id]

# Only interested in HLA-A, B, C
patient_alleles_dt <- patient_alleles_dt[gene %in% c("A", "B", "C")]

# check no genes are repeated
n_genes <- patient_alleles_dt[,.N,by = gene]
if(nrow(n_genes[N>1])>0){
  msg1 <- "These genes are repeated in the HLA alleles:\n"
  msg2 <- paste0(n_genes[N>1]$GENE, collapse = "\n")
  stop(msg1, msg2)
}

# check the patients HLAs are in the reference 
cat("Checking if HLA alleles are in the reference\n")

HLA_fasta <- read.fasta(mhc_fasta_path)
hla_fasta_alleles <- names(HLA_fasta)
patient_alleles <- unique(c(patient_alleles_dt$allele1, patient_alleles_dt$allele2))
if(any(!patient_alleles %in% hla_fasta_alleles)){
  stop("HLAHD has predicted alleles that are not in the reference - are you sure the input reference is from the same IMGT version as HLAHD used?")
}

##### Align HLA alleles #####

cat("Aligning HLA alleles\n")

patient_alleles_dt[, num_snps := as.numeric(NA)]
patient_alleles_dt[allele1 == allele2, homozygous := TRUE]
patient_alleles_dt[allele1 != allele2, homozygous := FALSE]

for(row_idx in 1:nrow(patient_alleles_dt)){
  
  if(patient_alleles_dt[row_idx]$homozygous){
    next
  }
  
  allele1 <- patient_alleles_dt[row_idx]$allele1
  allele2 <- patient_alleles_dt[row_idx]$allele2
  
  # the reference sequence for the patient's two HLA alleles
  allele1_fasta <- HLA_fasta[[allele1]]
  allele2_fasta <- HLA_fasta[[allele2]]
  
  cat("Aligning ", allele1, " and ", allele2, "\n")
  
  #perform local pairwise alignement
  allele1_seq <- paste0(toupper(as.character(allele1_fasta)), collapse ="")
  allele2_seq <- paste0(toupper(as.character(allele2_fasta)), collapse ="")
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  
  allele1_allele2_align <- Biostrings::pairwiseAlignment(allele1_seq, allele2_seq, 
                                                         substitutionMatrix = sigma,
                                                         gapOpening = -2,
                                                         gapExtension = -4,
                                                         scoreOnly = FALSE,
                                                         type='local')
  
  mismatch_positions_tmp <- getMisMatchPositionsPairwiseAlignment(allele1_allele2_align)
  
  # should be the same
  if(length(mismatch_positions_tmp$diffSeq1) != length(mismatch_positions_tmp$diffSeq2)){
    stop("What's going on?")
  }
  
  n_mismatch <- length(mismatch_positions_tmp$diffSeq1)
  
  patient_alleles_dt[row_idx, num_snps := n_mismatch]
  
  allele1_bed <- data.frame(Chrom = allele1,
                            chromStart = mismatch_positions_tmp$diffSeq1 - 1,
                            chromEnd = mismatch_positions_tmp$diffSeq1)
  
  allele1_bed_path <- paste0(allele1, "_", genome_transcriptome, ".snp_pos.bed")
  
  fwrite(allele1_bed, allele1_bed_path, sep = "\t", col.names = FALSE)
  
  
  allele2_bed <- data.frame(Chrom = allele2,
                            chromStart = mismatch_positions_tmp$diffSeq2 - 1,
                            chromEnd = mismatch_positions_tmp$diffSeq2)
  allele2_bed_path <- paste0(allele2, "_", genome_transcriptome, ".snp_pos.bed")
  
  fwrite(allele2_bed, allele2_bed_path, sep = "\t", col.names = FALSE)
}

patient_alleles_dt[homozygous == TRUE, num_snps := 0]

# write out the MHC_HAMMER_output table
fwrite(patient_alleles_dt, file = paste0(patient_id, "_", genome_transcriptome, "_allele_table.csv"))
