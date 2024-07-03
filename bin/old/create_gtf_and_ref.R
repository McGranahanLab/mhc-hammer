suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--gtf_out_path',  nargs=1,
                    help='Where to save gtf',
                    required=TRUE)
parser$add_argument('--fa_out_path',  nargs=1,
                    help='Where to save fa',
                    required=TRUE)
parser$add_argument('--hla_path',  nargs=1,
                    help="Path to individual's HLA alleles",
                    required=TRUE)
parser$add_argument('--mhc_fasta',  nargs=1,
                    help='Path to HLA fasta',
                    required=TRUE)
parser$add_argument('--mhc_gtf',  nargs=1,
                    help='Path to HLA gtf',
                    required=TRUE)
parser$add_argument('--genome_out_path',  nargs=1,
                    help='Path to genome_size table',
                    required=TRUE)

args <- parser$parse_args()
gtf_out_path <- args$gtf_out_path
fa_out_path <- args$fa_out_path
hla_path <- args$hla_path
hla_fasta_path <- args$mhc_fasta
genome_out_path <- args$genome_out_path
gtf_path <- args$mhc_gtf

hlahd_tab <- fread(hla_path, header = FALSE)
hla_alleles <- unique(c(hlahd_tab$V2, hlahd_tab$V3))

# remove any cases of "not typed"
hla_alleles <- hla_alleles[!grepl("not typed", hla_alleles)]

# get HLA fasta
HLA_fasta <- read.fasta(hla_fasta_path)
patient_fasta <- HLA_fasta[hla_alleles]

# get gtf
gtf <- fread(gtf_path)
patient_gtf <- gtf[seqname %in% hla_alleles]
# patient_gtf <- patient_gtf[type %in% c("5UTR", "3UTR", "exon")]
patient_gtf <- patient_gtf[order(seqname, feature_start)]

# save gtf file
fwrite(patient_gtf, 
       file = gtf_out_path,
       sep = "\t", col.names = FALSE, quote = FALSE)

# save fasta file
write.fasta(patient_fasta,
            file = fa_out_path,
            names = names(patient_fasta))

genome_length <- sum(sapply(patient_fasta, length))
genome_size <- data.frame(genome_size = min(c(floor(log(genome_length, base = 2)/2 - 1)), 14))
fwrite(genome_size, file = genome_out_path, col.names = FALSE, quote = FALSE, sep = "\t")
