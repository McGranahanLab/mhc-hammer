suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# path to rpkm tables - can be a list of paths
parser$add_argument('--rpkm_tables' , nargs='+',
                    help='Path to rpkm tables',
                    required=TRUE)

args <- parser$parse_args()

rpkm_tables_path <- args$rpkm_tables

# read in rpkm tables
rpkm_tables <- lapply(rpkm_tables_path, function(x) fread(x))

# rbind all rpkm tables
rpkm_tables <- rbindlist(rpkm_tables)

# make dt with just allele1 info
allele1_dt <- rpkm_tables[, .(patient, sample_name, gene, allele1, allele1_rpkm, 
                              allele1_updated_count,
                              homozygous, frac_map_uniquely, frac_reads_multi_gene)]

setnames(allele1_dt, 
         c("allele1", "allele1_rpkm", "allele1_updated_count"),
         c("allele", "allele_rpkm", "allele_updated_count"))

# do the same for allele2
allele2_dt <- rpkm_tables[, .(patient, sample_name, gene, allele2, allele2_rpkm, 
                              allele2_updated_count,
                              homozygous, frac_map_uniquely, frac_reads_multi_gene)]

setnames(allele2_dt, 
         c("allele2", "allele2_rpkm", "allele2_updated_count"),
         c("allele", "allele_rpkm", "allele_updated_count"))

# rbind the two dts
allele_rpkm_tables <- rbindlist(list(allele1_dt, allele2_dt))

# order by patient and sample_name and gene
allele_rpkm_tables <- allele_rpkm_tables[order(patient, sample_name, gene)]

# if a sample is homozygous, remove second allele
allele_rpkm_tables <- allele_rpkm_tables[homozygous == FALSE | (homozygous == TRUE & allele != "")]

# set column order
setcolorder(allele_rpkm_tables, 
            c("patient", "sample_name", "gene", "allele", "homozygous", 
              "frac_map_uniquely", "frac_reads_multi_gene", 
              "allele_updated_count", "allele_rpkm"))

# write out final allele rpkm table
fwrite(allele_rpkm_tables, file = "cohort_allele_rpkm.csv")

# Now make gene level rpkm table
gene_rpkm_table <- unique(rpkm_tables[, .(patient, sample_name, 
                                          allele1, allele2, gene, homozygous,
                                          frac_reads_multi_gene,
                                          gene_rpkm)])

# order by patient and sample_name and gene
gene_rpkm_table <- gene_rpkm_table[order(patient, sample_name, gene)]

# write out final gene rpkm table
fwrite(gene_rpkm_table, file = "cohort_gene_rpkm.csv")
