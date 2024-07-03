suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--hlahd_folder',  nargs=1,
                    help='Path to HLAHD output.',
                    required=TRUE)

parser$add_argument('--gtf_path',  nargs=1,
                    help='Path to GTF.',
                    required=TRUE)

parser$add_argument('--sample_id',  nargs=1,
                    help='sample_id',
                    required=TRUE)

parser$add_argument('--genes',  nargs='+',
                    help='genes to run',
                    required=TRUE)

args <- parser$parse_args()
hlahd_result_folder <- args$hlahd_folder
gtf_path <- args$gtf_path
sample_id <- args$sample_id
genes <- args$genes
hlahd_predictions_folder <- "hla_predictions"

# read the gtf and add columns that we'll use later
gtf <- fread(gtf_path)
gtf <- gtf[type %in% c("exon", "intron")]
gtf[, feature_seq_in_imgt := as.logical(gsub('.*?feature_seq_in_imgt "(.*?)";.*', "\\1", other_values))]
gtf[, similar_allele := gsub('.*?similar_allele_save_name "(.*?)";.*', "\\1", other_values)]
gtf[, similar_allele_distance := gsub('.*?min_dist "(.*?)";.*', "\\1", other_values)]

# get the HLA-HD files and turn into one table
est_files <- paste0(hlahd_result_folder, "/", sample_id, "_", genes, ".est.txt")
est_files_exist <- file.exists(est_files)

if(any(!est_files_exist)){
  cat("These HLA-HD *est.txt files are missing:\n")
  cat(est_files[!est_files_exist], sep = "\n")
  stop()
}

hlahd_result <- data.table(gene = character(),
                     allele_id = character(),
                     allele = character(),
                     digit2 = character(),
                     digit4 = character(),
                     digit6 = character(),
                     digit8 = character())

for(g_file in est_files){
  cat("doing: ", g_file, "\n")
  g <- readLines(g_file)
  
  if(g[1] == "No candidate."){
    next()
  }
  
  alleles <- g[3]
  alleles <- strsplit(x = alleles, split = "\t")[[1]]
  
  # allele1
  allele1 <- alleles[1]
  allele1 <- strsplit(allele1, ",")[[1]]
  
  allele1_fields <- gsub(".*\\*", "", allele1)
  allele1_fields <- strsplit(allele1_fields, ":")
  allele1_fields_max_length <- max(sapply(allele1_fields, FUN = length))
  allele1_fields <- lapply(allele1_fields, function(x) c(x, rep("NA", allele1_fields_max_length - length(x))) )
  
  allele1_fields_dt <- transpose(as.data.table(allele1_fields))
  setnames(allele1_fields_dt, 
           old = paste0("V", 1:ncol(allele1_fields_dt)),
           new = paste0("digit", seq(from = 2, to = ncol(allele1_fields_dt)*2, by = 2)))
  
  allele1_fields_dt[, gene := gsub(".*_(.*).est.txt", "\\1", basename(g_file))]
  allele1_fields_dt[, allele_id := "allele1"]        
  allele1_fields_dt[, allele := allele1]        
  
  hlahd_result <- rbindlist(list(hlahd_result, allele1_fields_dt), use.names = TRUE, fill = TRUE)
  
  # allele2
  allele2 <- alleles[2]
  if(allele2 != "-"){
    allele2 <- strsplit(allele2, ",")[[1]]
    
    allele2_fields <- gsub(".*\\*", "", allele2)
    allele2_fields <- strsplit(allele2_fields, ":")
    allele2_fields_max_length <- max(sapply(allele2_fields, FUN = length))
    allele2_fields <- lapply(allele2_fields, function(x) c(x, rep("NA", allele2_fields_max_length - length(x))) )
    
    allele2_fields_dt <- transpose(as.data.table(allele2_fields))
    setnames(allele2_fields_dt, 
             old = paste0("V", 1:ncol(allele2_fields_dt)),
             new = paste0("digit", seq(from = 2, to = ncol(allele2_fields_dt)*2, by = 2)))
    
    allele2_fields_dt[, gene := gsub(".*_(.*).est.txt", "\\1", basename(g_file))]
    allele2_fields_dt[, allele_id := "allele2"]        
    allele2_fields_dt[, allele := allele2]       
    hlahd_result <- rbindlist(list(hlahd_result, allele2_fields_dt), use.names = TRUE, fill = TRUE)
  }else{
    allele1_fields_dt[,allele_id := "allele2"]
    hlahd_result <- rbindlist(list(hlahd_result, allele1_fields_dt), use.names = TRUE, fill = TRUE)
  }
}

# Check there is a HLAHD result for at least one gene
any_gene_present = FALSE

for (gene_id in genes) {
  gene_data = hlahd_result[gene == gene_id]
  
  if (all(c("allele1", "allele2") %in% gene_data$allele_id)) {
    any_gene_present = TRUE
    break
  }
}

# Create HLAHD output when there is at least one gene present
if(any_gene_present){

  # add in the number of exons with sequence and the number of introns with a sequence in IMGT
  n_exons_with_seq_dt <- gtf[type == "exon", sum(feature_seq_in_imgt), by = "seqname"]
  setnames(n_exons_with_seq_dt, c("seqname", "n_exons_with_seq"))
  n_exons_dt <- gtf[type == "exon", .N, by = "seqname"]
  setnames(n_exons_dt, c("seqname", "n_exons"))

  n_introns_with_seq_dt <- gtf[type == "intron", sum(feature_seq_in_imgt), by = "seqname"]
  setnames(n_introns_with_seq_dt, c("seqname", "n_introns_with_seq"))
  n_introns_dt <- gtf[type == "intron", .N, by = "seqname"]
  setnames(n_introns_dt, c("seqname", "n_introns"))

  exon_introns_with_seq_dt <- merge(n_exons_with_seq_dt, n_exons_dt, by = "seqname", all = TRUE)
  exon_introns_with_seq_dt <- merge(exon_introns_with_seq_dt, n_introns_with_seq_dt, by = "seqname", all = TRUE)
  exon_introns_with_seq_dt <- merge(exon_introns_with_seq_dt, n_introns_dt, by = "seqname", all = TRUE)
  exon_introns_with_seq_dt[,exon_frac := paste0(n_exons_with_seq, "/", n_exons)]
  exon_introns_with_seq_dt[,intron_frac := paste0(n_introns_with_seq, "/", n_introns)]

  # add in if there is a full transcriptome seq and genome seq
  exon_introns_with_seq_dt[n_exons_with_seq == n_exons & n_introns_with_seq == n_introns, full_genome_seq_in_imgt := TRUE]
  exon_introns_with_seq_dt[is.na(full_genome_seq_in_imgt), full_genome_seq_in_imgt := FALSE]

  exon_introns_with_seq_dt[n_exons_with_seq == n_exons, full_transcriptome_seq_in_imgt := TRUE]
  exon_introns_with_seq_dt[is.na(full_transcriptome_seq_in_imgt), full_transcriptome_seq_in_imgt := FALSE]

  # add the most similar allele and the distance
  similar_allele <- unique(gtf[,c("seqname", "similar_allele", "similar_allele_distance")])
  exon_introns_with_seq_dt <- merge(exon_introns_with_seq_dt, similar_allele, by = "seqname", all = TRUE)
  exon_introns_with_seq_dt[similar_allele == seqname, similar_allele := NA]

  # add this information to the hlahd results
  hlahd_result[, seqname := tolower(gsub(":|\\*|-", "_", allele))]
  hlahd_result <- merge(hlahd_result, exon_introns_with_seq_dt, by = "seqname", all.x = TRUE)

  if(nrow(hlahd_result[is.na(full_transcriptome_seq_in_imgt)]) > 0){
    stop("Why are there NAs in this table? Does the GTF come from the same IMGT version that HLA-HD used?")
  }

  hlahd_gene_estimates_path <- paste0(hlahd_result_folder, "/hlahd_gene_estimates.csv")
  fwrite(hlahd_result, hlahd_gene_estimates_path)

  # order alleles by how many exons and introns have sequence in IMGT
  hlahd_result <- hlahd_result[order(-n_exons_with_seq, -n_introns_with_seq)]

  # choose the top alleles for each gene
  hlahd_result[,gene_idx := 1:.N, by = c("gene", "allele_id")]
  hlahd_result <- hlahd_result[gene_idx == 1]

  # make wide table
  hlahd_result_wide <- dcast(hlahd_result, gene ~ allele_id, value.var = "seqname")

  # If genes are missing in the gene columm, add a row with "not typed" in the allele columns
  if(!all(genes %in% hlahd_result_wide$gene)){
    missing_genes <- genes[!genes %in% hlahd_result_wide$gene]
    missing_genes <- data.table(gene = missing_genes, allele1 = "not typed", allele2 = "not typed")
    hlahd_result_wide <- rbind(hlahd_result_wide, missing_genes)
    # arrange alphabetically by gene
    hlahd_result_wide <- hlahd_result_wide[order(gene)]
  }

  hla_alleles_path <- paste0(hlahd_result_folder, "/", sample_id, "_hla_alleles.csv")
  cat("Saving to: ", hla_alleles_path, "\n")
  fwrite(hlahd_result_wide, hla_alleles_path, col.names = FALSE)  

} else { # if no HLA genes were found, create a file with "not typed" in the allele columns
  hla_alleles_path <- paste0(hlahd_result_folder, "/", sample_id, "_hla_alleles.csv")
  cat("Saving to: ", hla_alleles_path, "\n")
  fwrite(data.table(gene = genes, allele1 = "not typed", allele2 = "not typed"), hla_alleles_path, col.names = FALSE)
}