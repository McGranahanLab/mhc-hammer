library(data.table)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('--gtf_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--gtf_vep_path',  nargs=1,
                    required=TRUE)
parser$add_argument('--allele',  nargs=1,
                    required=TRUE)

args <- parser$parse_args()
gtf_path <- args$gtf_path
# gtf_path <- "/nemo/project/proj-tracerx-lung/tctProjects/putticc/mhc_pipeline_mutations/LOHHLA/test_muts/input_data/test.gtf"
gtf_vep_path <- args$gtf_vep_path
allele <- args$allele
# allele <- "hla_b_57_01_01_01"
gene_name <- toupper(gsub("(hla_.*?)_.*", "\\1", allele))
gene_name <- gsub("_", "-", gene_name)

gtf <- fread(gtf_path)
gtf <- gtf[V1 == allele]

# make exon table
exon_dt <- gtf[V3 == "exon"]
exon_dt[,exon_number := gsub('.*exon_number "(.*?)".*', "\\1", V9)]
exon_dt[,exon_name := gsub('.*feature_name "(.*?)".*', "\\1", V9)]
exon_dt[,new_final_col := paste0('gene_id "', gene_name, '"; transcript_id "', allele, 
                                 '"; exon_number "exon_', exon_number, '"; exon_id "',exon_name, '";')]
exon_dt <- exon_dt[,c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "new_final_col")]

# add cds
cds_dt <- copy(exon_dt)
cds_dt[,V3 := "CDS"]
cds_dt[,new_final_col := gsub("exon_id", "ccds_id", new_final_col)]


# add gene and transcript line
gene_dt <- data.table(V1 = allele,
                      V2 = "imgt",
                      V3 = "gene",
                      V4 = min(gtf$V4),
                      V5 = max(gtf$V5),
                      V6 = unique(gtf$V6),
                      V7 = unique(gtf$V7),
                      V8 = unique(gtf$V8),
                      new_final_col = paste0('gene_id "', gene_name, '"; gene_name "', gene_name, '";'))

transcript_dt <- data.table(V1 = allele,
                      V2 = "imgt",
                      V3 = "transcript",
                      V4 = min(exon_dt$V4),
                      V5 = max(exon_dt$V5),
                      V6 = unique(gtf$V6),
                      V7 = unique(gtf$V7),
                      V8 = unique(gtf$V8),
                      new_final_col = paste0('gene_id "', gene_name, '"; transcript_id "', allele, '"; gene_name "',
                                  gene_name, '"; transcript_name "', allele, '"; transcript_biotype "protein_coding";'))

# put together and sort
vep_gtf <- rbindlist(list(gene_dt, transcript_dt, exon_dt, cds_dt), use.names = TRUE)
vep_gtf <- vep_gtf[order(V4)]

fwrite(vep_gtf, file = gtf_vep_path, sep = "\t", col.names = FALSE, quote = FALSE)

