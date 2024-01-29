library(seqinr)
library(data.table)
library(argparse)

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--imgt_dir', nargs=1,
                    help='Path to IMGT directory.',
                    required=TRUE)
parser$add_argument('--out_dir', nargs=1,
                    help='Directory to save output.',
                    required=TRUE)

args <- parser$parse_args()
imgt_dir <- args$imgt_dir
out_dir <- args$out_dir
reference_type <- args$reference_type

# make output dir 
cat('Creating output directory: ', out_dir, '\n')
if(!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}else{
  cat('Directory already exists.\n')
}

# check IMGT dir exists
if(!dir.exists(imgt_dir)){
  stop("IMGT directory does not exist: ", imgt_dir)
}

fasta_dir <- paste0(imgt_dir, "/fasta")
fasta_files <- list.files(fasta_dir)

gen_fasta_files <- grep("gen", fasta_files, value = TRUE)
nuc_fasta_files <- grep("nuc", fasta_files, value = TRUE)

cat('Reading fasta files\n')
all_fasta <- c()
for(f in c(gen_fasta_files, nuc_fasta_files)){
  f_path <- paste0(fasta_dir, "/", f)
  cat(f_path, "\n")
  all_fasta <- c(all_fasta, read.fasta(f_path, forceDNAtolower = FALSE))
}

rev_all_fasta <- vector(mode = 'list', length = length(all_fasta))
for(list_idx in 1:length(all_fasta)){
 
  seq <- as.character(all_fasta[[list_idx]])
  
  g_idx <- which(seq == "G")
  c_idx <- which(seq == "C")
  t_idx <- which(seq == "T")
  a_idx <- which(seq == "A")
  
  seq[g_idx] <- "C"
  seq[c_idx] <- "G"
  seq[t_idx] <- "A"
  seq[a_idx] <- "T"
  
  rev_all_fasta[[list_idx]] <- rev(seq)
  
}

names(rev_all_fasta) <- paste0(names(all_fasta), "_rev")

forward_rev_fasta <- c(rev_all_fasta, all_fasta)

cat('Saving fasta file\n')
all_fasta_path <- paste0(out_dir, "/all_fasta.fasta")
write.fasta(forward_rev_fasta,
            file=paste0(out_dir, "/all_fasta.fasta"),
            names = names(forward_rev_fasta),nbchar=80)

