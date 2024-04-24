library(data.table)
library(seqinr)
library(argparse)
# https://github.com/ANHIG/IMGTHLA/blob/Latest/Manual.md
parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--path_to_hla_dat', nargs=1,
                    help='Path to IMGT hla.dat',
                    required=TRUE)
parser$add_argument('--save_dir', nargs=1,
                    help='Path directory to save the output files in',
                    required=TRUE)
parser$add_argument('--functions_file', nargs=1,
                    help='Path to the file mhc_reference_functions.R',
                    required=TRUE)

args <- parser$parse_args()
path_to_hla_dat <- args$path_to_hla_dat
save_dir <- args$save_dir
functions_file <- args$functions_file

source(functions_file)

##### process hla dat #####

cat("Reading in the hla.dat file\n")
x <- readLines(path_to_hla_dat)

# get the id lines
id_lines <- grep('^ID', x)
end_lines <- id_lines[2:length(id_lines)] -1
end_lines <- c(end_lines, length(x))
allele_index <- data.table(start_idx = id_lines,
                           end_idx = end_lines)

all_allele_features <- data.table()
all_allele_info <- data.table()
alleles_with_missing_data <- data.table()

cat("Processing the hla.dat file\n")
for(line_idx in 1:nrow(allele_index)){
  
  if( (line_idx %% 1000) == 0){
    cat("Processed: ", line_idx, "/", nrow(allele_index), "alleles\n")
  }else if( line_idx == nrow(allele_index) ){
    cat("Processed: ", line_idx, "/", nrow(allele_index), "alleles\n")
  }
  
  allele_entry <- x[allele_index[line_idx]$start_idx : allele_index[line_idx]$end_idx]
  
  id_line <- grep('^ID', allele_entry, value = TRUE)
  de_line <- grep('^DE', allele_entry, value = TRUE)
  kw_line <- grep('^KW', allele_entry, value = TRUE)
  ft_lines <- grep('^FT', allele_entry, value = TRUE)
  sq_line <- grep('^SQ', allele_entry, value = TRUE)
  
  if(length(ft_lines) == 0){
    cat("No FT line for line_idx: ",  line_idx, "\n")
    next
  }
  
  # get cell line data if it exists
  cell_line_lines <- grep("cell_line", ft_lines, value = TRUE)
  cell_lines <- paste0(gsub('.*cell_line=\\\"(.*)\\\"', "\\1", cell_line_lines), collapse = ",")
  
  # get codon start data if it exists
  codon_start_lines <- grep("codon_start=", ft_lines, value = TRUE)
  codon_start <- gsub('.*codon_start=(.*)', "\\1", codon_start_lines)
  
  # get gene
  gene_line <- grep("gene=", ft_lines, value = TRUE)
  gene <- paste0(gsub('.*gene=\\\"(.*)\\\"', "\\1", gene_line), collapse = ",")
  
  if(!gene %in% c("HLA-A", "HLA-B", "HLA-C")){
    next
  }
  
  # get allele
  allele_line <- grep("allele=", ft_lines, value = TRUE)
  allele <- paste0(gsub('.*allele=\\\"(.*)\\\"', "\\1", allele_line), collapse = ",")
  
  # get ethnic
  ethnic_line <- grep("ethnic=", ft_lines, value = TRUE)
  ethnic <- paste0(gsub('.*ethnic=\\\"(.*)\\\"', "\\1", ethnic_line), collapse = ",")
  
  # check if it is partial
  partial <- any(grepl("FT[[:space:]].*/partial", ft_lines))
  
  # get mol type
  mol_type_line <- grep("mol_type", ft_lines, value = TRUE)
  mol_type <- gsub('.*mol_type=\"(.*)\"', "\\1", mol_type_line)
  
  cds_line <- get_cds_line(ft_lines)
  
  # get the translation line
  translation <- get_translation_line(ft_lines)
  
  # get coordinates line
  coords_dt <- get_coords(ft_lines)
  
  # get sequence length
  seq_length <- as.numeric(gsub(".* Sequence ([0-9]+).*", "\\1", sq_line))
  
  # get sequence
  sq_line_idx <- grep('^SQ', allele_entry)
  seq <- c()
  for(i in (sq_line_idx + 1) : length(allele_entry)){
    if(allele_entry[i] =="//" ){
      break
    }
    seq_line <-  gsub("([a-z]+.*?)[0-9].*", "\\1", allele_entry[i])
    seq <- c(seq, gsub("[[:space:]]", "", seq_line))
  }
  seq <- paste0(seq, collapse = "")
  if(nchar(seq) != seq_length){
    stop("Sequence length and the number of bases is not the same")
  }
  
  if(nchar(seq) != seq_length){
    stop("Why do the seq lengths differ?")
  }
  
  all_allele_info <- rbindlist(list(all_allele_info,
                                    data.table(allele_name = allele,
                                               cell_lines = cell_lines,
                                               gene = gene,
                                               ethnic = ethnic,
                                               mol_type = mol_type,
                                               partial = partial,
                                               codon_start = codon_start,
                                               cds_line = cds_line,
                                               translation = translation,
                                               seq_length = seq_length,
                                               seq = toupper(seq))))
  
  # add the feature sequence to the coords dt
  coords_dt[,length := as.numeric(end) - as.numeric(start) + 1]
  coords_dt[,allele_name := allele]
  coords_dt[,gene := gene]
  if(sum(coords_dt$length) != seq_length){
    stop("Sequence length and the seq coords is not the same")
  }
  
  for(i in 1:nrow(coords_dt)){
    coords_dt[i,feature_seq :=  substr(seq, start = coords_dt[i]$start, stop = coords_dt[i]$end)]
  }
  
  # add in the CDS sequence
  # if(!grepl("join", cds_line)){
  #   alleles_with_missing_data <- rbindlist(list(alleles_with_missing_data,
  #                                               data.table(allele_name = allele,
  #                                                          type = "cds_no_join",
  #                                                          line_idx = line_idx)))
  #   next
  # }
  
  cds_coords <- gsub(".*?([0-9].*[0-9]).*$", "\\1", cds_line)
  # cds_coords <- gsub('join\\(', "", cds_line)
  # cds_coords <- gsub(')', "", cds_coords)
  cds_dt <- data.table(start_end = strsplit(cds_coords, ",")[[1]])
  cds_dt[,start := gsub("\\.\\..*", "", start_end)]
  cds_dt[,end := gsub(".*\\.\\.", "", start_end)]
  cds_dt[,start_end := NULL]
  cds_dt[, start := suppressWarnings(as.numeric(start))]
  cds_dt[, end := suppressWarnings(as.numeric(end))]
  
  if( nrow(cds_dt[is.na(start)]) > 0 | nrow(cds_dt[is.na(end)]) > 0 ){
    
    stop("start or stop of CDS not numeric?\n")
    # coords_dt[,feature_seq := toupper(feature_seq)]
    # all_allele_features <- rbindlist(list(all_allele_features, coords_dt), use.names = TRUE)
    # alleles_with_missing_data <- rbindlist(list(alleles_with_missing_data,
    #                                             data.table(allele_name = allele,
    #                                                        type = "cds_start_stop_not_numeric",
    #                                                        line_idx = line_idx)))
    
  }
  
  # double check that the CDS isnt more than the length of the sequence
  if(max(c(cds_dt$end, cds_dt$start)) > seq_length){
    stop("CDS coordinates higher than sequence legnth")
  }
  
  for(i in 1:nrow(cds_dt)){
    cds_dt[i,feature_seq :=  substr(seq, start = cds_dt[i]$start, stop = cds_dt[i]$end)]
  }
  
  cds_dt[,length := as.numeric(end) - as.numeric(start) + 1]
  cds_dt[,allele_name := allele]
  cds_dt[,gene := gene]
  cds_dt[,type := "CDS"]
  
  coords_dt <- rbindlist(list(coords_dt, cds_dt), use.names = TRUE, fill = TRUE)
  coords_dt[,feature_seq := toupper(feature_seq)]
  all_allele_features <- rbindlist(list(all_allele_features, coords_dt), use.names = TRUE)
  
  
}

all_allele_features[, allele_name := gsub("HLA-", "", allele_name)]
all_allele_info[, allele_name := gsub("HLA-", "", allele_name)]

# make output dir 
cat('Creating output directory: ', save_dir, '\n')
if(!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
}else{
  cat('Directory already exists.\n')
}

cat("Saving output files\n")
fwrite(all_allele_features, 
       file = paste0(save_dir, "/all_allele_features.csv"))
fwrite(all_allele_info, 
       file = paste0(save_dir, "/all_allele_info.csv"))
