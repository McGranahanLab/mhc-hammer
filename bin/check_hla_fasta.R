library(data.table)
library(seqinr)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('--hla_fasta_path',  nargs=1,
                  help='Path to HLA fasta file.',
                  required=TRUE)

parser$add_argument('--hla_allele_files', nargs = '+',
                  help="Paths to HLA allele files",
                  required=TRUE)

parser$add_argument('--parsed_csv', nargs = 1,
                  help="parsed csv file",
                  required=TRUE)

args <- parser$parse_args()

# Read in the fasta file
fasta_path <- args$hla_fasta_path
hla_fasta <- read.fasta(fasta_path)

# Read in the parsed csv file
parsed_csv_path <- args$parsed_csv
parsed_csv <- fread(parsed_csv_path)

# Read in the allele files
hla_allele_paths <- args$hla_allele_files
hla_alleles_long <- lapply(seq_along(hla_allele_paths), function(i) {
  dt <- fread(hla_allele_paths[i], sep = ",", header = FALSE, col.names = c("gene", "allele_1", "allele_2"))
  file_name <- basename(hla_allele_paths[i])  # Extract file name
  melt(dt, id.vars = "gene", measure.vars = c("allele_1", "allele_2"), 
       value.name = "allele")[, file := file_name][]  # Reshape and add file column
})

combined_hla_alleles_long <- rbindlist(hla_alleles_long)

# check which alleles are missing from the fasta file
missing_alleles <- combined_hla_alleles_long[!allele %in% names(hla_fasta),]

# if there are any missing alleles, stop and print the missing alleles
if(nrow(missing_alleles) > 0){
  missing_alleles <- missing_alleles[,.(gene, allele, file)]
  stop("The custom HLA alleles provided are not present within the current HLA allele list\nAre you using the correct IMGT version / the correct hla allele formatting? Check the README for more details!\nThe following alleles are missing from the HLA fasta file:\n", paste0(missing_alleles$allele, " (", missing_alleles$file, ")", collapse = "\n"), "\n")
}

# write out parsed csv as it is now validated
fwrite(parsed_csv, file = "samplesheet.valid.csv")

