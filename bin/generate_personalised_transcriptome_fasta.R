suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--patient_id',  nargs=1,
                    help='Patient ID.',
                    required=TRUE)

parser$add_argument('--patient_hla_alleles',  nargs=1,
                    help='Patient HLA alleles file.',
                    required=TRUE)

parser$add_argument('--mhc_fasta_path',  nargs=1,
                    help='HLA fasta path.',
                    required=TRUE)

args <- parser$parse_args()

patient_id <- args$patient_id
patient_alleles <- args$patient_hla_alleles
HLA_fasta_path <- args$mhc_fasta_path

msg <- paste0("generating personalised transcriptome fasta")
print(msg)

# read in hla fasta
HLA_fasta <- read.fasta(HLA_fasta_path)

# read in patient hla alleles table
hla_alleles_dt <- fread(patient_alleles)
hla_alleles <- unique(c(hla_alleles_dt$allele1, hla_alleles_dt$allele2))

# remove any cases of "not typed"
hla_alleles <- hla_alleles[!grepl("not typed", hla_alleles)]

# get patient hla fasta 
patient.HLA_fasta <- HLA_fasta[hla_alleles]

# write patient hla fasta 
patient.hlaFasta_path <- paste0(patient_id, '_mhc_transcriptome_reference.fa')
write.fasta(patient.HLA_fasta,
            file = patient.hlaFasta_path,
            names = names(patient.HLA_fasta))

print("Done!")
