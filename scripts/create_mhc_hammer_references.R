library(data.table)
library(seqinr)
library(argparse)

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--imgt_dir', nargs=1,
                    help='Path to IMGT dir.',
                    required=TRUE)

parser$add_argument('--all_allele_features_path', nargs=1,
                    help='Path to hla_dat.',
                    required=TRUE)

parser$add_argument('--strand_info', nargs=1,
                    required=TRUE)

parser$add_argument('--outdir', nargs=1,
                    help='Path to hla_dat.',
                    required=TRUE)

parser$add_argument('--functions_file', nargs=1,
                    help='Path to the file mhc_reference_functions.R',
                    required=TRUE)

args <- parser$parse_args()
imgt_dir <- args$imgt_dir
all_allele_features_path <- args$all_allele_features_path
strand_info <- args$strand_info
outdir <- args$outdir
functions_file <- args$functions_file

source(functions_file)

cat("Getting genome fasta\n")
gen_fasta_path <- paste0(imgt_dir, "/hla_gen.fasta")
gen_fasta <- read.fasta(gen_fasta_path, forceDNAtolower = FALSE)
allele_names <- NULL
for(i in 1:length(gen_fasta)){allele_names = c(allele_names, attr(gen_fasta[[i]], "Annot"))}
allele_names <- unname(sapply(allele_names, function(x) unlist(strsplit(x, " "))[2]))
names(gen_fasta) <- allele_names

cat("Getting cds fasta\n")
nuc_fasta_path <- paste0(imgt_dir, "/hla_nuc.fasta")
nuc_fasta <- read.fasta(nuc_fasta_path, forceDNAtolower = FALSE)
allele_names <- NULL
for(i in 1:length(nuc_fasta)){allele_names = c(allele_names, attr(nuc_fasta[[i]], "Annot"))}
allele_names <- unname(sapply(allele_names, function(x) unlist(strsplit(x, " "))[2]))
names(nuc_fasta) <- allele_names

cat("Getting allele features table\n")
hla_dat_allele_features <- fread(all_allele_features_path)

all_allele_list <- get_allele_list(paste0(imgt_dir, "/Allele_status.txt"))

cat("Getting strand information\n")
strand_dt <- fread(strand_info, header = TRUE)

genome_fasta <- NULL
genome_strand_specific_fasta <- NULL
transcriptome_fasta <- NULL
gtf <- data.table()
exon_intron_dt <- data.table()
alleles_with_unmatching_cds_exon_seq <- c()
alleles_with_unmatching_gen_seq <- c() 
alleles_with_unmatching_cds_seq <- c()

for(gene_to_run in c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G")){
# for(gene_to_run in unique(all_allele_list$gene)){
  
  cat("Processing ", gene_to_run, "\n")
  
  allele_list <- all_allele_list[gene == gene_to_run]
  if(nrow(allele_list) == 0){
    stop("Gene not in allele list ")
  }
  
  strand <- as.character(strand_dt[gene == gene_to_run]$strand)
  
  cat("\tGetting msf\n")
  msf_dir <- paste0(imgt_dir, "/msf")
  # get the distance matrix
  dist_mat <- get_dist_matrix(msf_dir, gene_to_run)
  # Update the distance matrices so it only includes alleles with full gdna sequences
  dist_mat <- dist_mat[,allele_list[Partial == "Full" & Type == "gDNA"]$Allele]
  
  cat("\tCreating references\n")
  
  for(line_idx in 1:nrow(allele_list)){
    
    if( (line_idx %% 1000) == 0){
      cat("\tProcessed: ", line_idx, "/", nrow(allele_list), "alleles\n")
    }else if( line_idx == nrow(allele_list) ){
      cat("\tProcessed: ", line_idx, "/", nrow(allele_list), "alleles\n")
    }
    
    allele <- allele_list[line_idx]$Allele
    save_allele_name <- paste0("hla_", (tolower(gsub('\\*|:', '_', allele))))
    
    ##### get genomic sequence ####
    if(allele_list[line_idx]$full_gen_seq_exists){
      
      # full genome sequence already exists
      allele_genome_seq <- as.character(gen_fasta[allele][[1]])
      if(length(allele_genome_seq) == 0){
        stop("Allele ", allele, " not in genome fasta")
      }
      allele_gen_fasta = as.SeqFastadna(allele_genome_seq)
      attributes(allele_gen_fasta)$name <- save_allele_name
      
      # we need this to make the gtf later
      genome_allele_dt <- hla_dat_allele_features[allele_name == allele & type != "CDS"]
      genome_allele_dt[,feature_seq_in_imgt := TRUE]
      genome_allele_dt[, min_dist := 0]
      genome_allele_dt[, similar_allele := allele]
      
      # double check it is the same as what is in the hla_dat file
      if(toupper(paste0(genome_allele_dt$feature_seq, collapse = "")) != paste0(allele_genome_seq, collapse = "")){
        # stop("The genome sequence from the fasta and the hla_dat aren't the same")
        alleles_with_unmatching_gen_seq <- c(alleles_with_unmatching_gen_seq, allele)
      }
      
    }else{
      # full genome sequence doesn't exist
      
      # get similar allele
      min_dist <- min(dist_mat[allele,])
      similar_allele <- names(which.min(dist_mat[allele,]))
      
      # get missing sequences
      genome_allele_dt <- get_missing_seq(allele, similar_allele, hla_dat_allele_features)
      genome_allele_dt[, min_dist := min_dist]
      genome_allele_dt[, similar_allele := similar_allele]
      
      # get the gen sequence
      allele_genome_seq <- unlist(strsplit(genome_allele_dt$feature_seq, split = ""))
      allele_gen_fasta = as.SeqFastadna(allele_genome_seq)
      attributes(allele_gen_fasta)$name <- save_allele_name
      
    }
    
    ##### get nuc sequence #####
    if(allele_list[line_idx]$full_nuc_seq_exists){
      
      # full nuclear sequence exists - this will be the CDS sequence
      allele_nuc_seq <- as.character(nuc_fasta[allele][[1]])
      if(length(allele_nuc_seq) == 0){
        stop("Allele ", allele, " not in nuc fasta")
      }
      
      # check whether the cds and the exon sequence is the same
      # if different we use the exon sequence
      allele_exon_seq_from_hla_dat <- toupper(paste0(hla_dat_allele_features[allele_name == allele & type == "exon"]$feature_seq, collapse = ""))
      allele_cds_seq_from_hla_dat <- toupper(paste0(hla_dat_allele_features[allele_name == allele & type == "CDS"]$feature_seq, collapse = ""))
      
      if(allele_cds_seq_from_hla_dat != paste0(allele_nuc_seq, collapse = "")){
        # cat("The CDS sequence and the exon sequence aren't the same for allele ", allele, ", will use the exon sequence for the transcriptome reference\n")
        alleles_with_unmatching_cds_seq <- c(alleles_with_unmatching_cds_seq, allele)
      }
      
      if(allele_exon_seq_from_hla_dat != paste0(allele_nuc_seq, collapse = "")){
        # cat("The CDS sequence and the exon sequence aren't the same for allele ", allele, ", will use the exon sequence for the transcriptome reference\n")
        alleles_with_unmatching_cds_exon_seq <- c(alleles_with_unmatching_cds_exon_seq, allele)
        # if different we use the exon sequence
        allele_nuc_seq <- allele_exon_seq_from_hla_dat
      }
      
      allele_nuc_fasta = as.SeqFastadna(allele_nuc_seq)
      attributes(allele_nuc_fasta)$name <- save_allele_name
      
    }else{
      
      # full nuc sequence doesnt exist
      
      # get similar allele
      min_dist <- min(dist_mat[allele,])
      similar_allele <- names(which.min(dist_mat[allele,]))
      
      # get missing sequences
      nuc_allele_dt <- get_missing_seq(allele, similar_allele, hla_dat_allele_features)
      
      # get the nuc sequence
      allele_nuc_seq <- unlist(strsplit(nuc_allele_dt[type == "exon"]$feature_seq, split = ""))
      allele_nuc_fasta = as.SeqFastadna(allele_nuc_seq)
      attributes(allele_nuc_fasta)$name <- save_allele_name
    }
    
    ##### get strand specific fasta #####
    # if the gene is on the reverse strand, this will reverse the sequence and change to the complement
    allele_gen_strand_specific_seq <- get_gen_seq_strand_specific(allele_genome_seq, strand)
    allele_gen_strand_specific_fasta = as.SeqFastadna(allele_gen_strand_specific_seq)
    attributes(allele_gen_strand_specific_fasta)$name <- save_allele_name
    
    ##### get gtf #####
    allele_gtf <- get_gtf(genome_allele_dt, save_allele_name, gene_to_run, strand)
    
    ##### combine the tables #####
    gtf <- rbindlist(list(gtf, allele_gtf), use.names = TRUE)
    genome_fasta[[save_allele_name]] <- allele_gen_fasta
    genome_strand_specific_fasta[[save_allele_name]] <- allele_gen_strand_specific_fasta
    transcriptome_fasta[[save_allele_name]] <- allele_nuc_fasta
  }
}

cat("Saving files")

if(!dir.exists(paste0(outdir, "/gtf/"))){
  dir.create(paste0(outdir, "/gtf/"), recursive = TRUE)
}

gtf_path <- paste0(outdir, "/gtf/mhc.gtf")
fwrite(gtf, file = gtf_path, sep = "\t", quote = FALSE)

if(!dir.exists(paste0(outdir, "/genome/"))){
  dir.create(paste0(outdir, "/genome/"), recursive = TRUE)
}

genome_path <- paste0(outdir, "/genome/mhc_genome.fasta")

write.fasta(genome_fasta,
            file = genome_path,
            names = names(genome_fasta), nbchar = 80)

genome_strand_specific_path <- paste0(outdir, "/genome/mhc_genome_strand.fasta")
write.fasta(genome_strand_specific_fasta,
            file = genome_strand_specific_path,
            names = names(genome_strand_specific_fasta), nbchar = 80)

if(!dir.exists(paste0(outdir, "/transcriptome/"))){
  dir.create(paste0(outdir, "/transcriptome/"), recursive = TRUE)
}

transcriptome_path <- paste0(outdir, "/transcriptome/mhc_cds.fasta")
write.fasta(transcriptome_fasta,
            file = transcriptome_path,
            names = names(transcriptome_fasta), nbchar = 80)

# save alleles_with_unmatching_nuc_seq
fwrite(data.table(allele = alleles_with_unmatching_cds_exon_seq),
       file = paste0(outdir, "/transcriptome/alleles_with_unmatching_cds_exon_seq.csv"))
fwrite(data.table(allele = alleles_with_unmatching_cds_seq),
       file = paste0(outdir, "/transcriptome/alleles_with_unmatching_cds_seq.csv"))
fwrite(data.table(allele = alleles_with_unmatching_gen_seq),
       file = paste0(outdir, "/transcriptome/alleles_with_unmatching_gen_seq.csv"))
