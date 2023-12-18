suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--hla_bam_paths',
                    help='Path to all HLA bams for a given sample.',
                    nargs="+",
                    required=TRUE)
parser$add_argument('--sample_name',
                    help='Sample name.',
                    required=TRUE)
parser$add_argument('--patient',
                    help='Sample name.',
                    required=TRUE)
parser$add_argument('--reference_fasta_path',
                    help='Path to reference.',
                    required=TRUE)
parser$add_argument('--out_file_prefix',
                    help='Prefix to save file as',
                    required=TRUE)

args <- parser$parse_args()
hla_bam_paths <- args$hla_bam_paths
sample_name <- args$sample_name
patient <- args$patient
reference_fasta_path <- args$reference_fasta_path
out_file_prefix <- args$out_file_prefix

hla_bam_paths <- grep("bam$", hla_bam_paths, value = TRUE)

make_read_table <- function(bam_path){
  
  bam <- scanBam(bam_path)
  
  reads_dt <- data.table(read = bam[[1]]$qname,
                         flag = bam[[1]]$flag,
                         chrom = bam[[1]]$rname,
                         mate_chrom = bam[[1]]$mrnm,
                         pos = bam[[1]]$pos,
                         strand = bam[[1]]$strand,
                         qwidth = bam[[1]]$qwidth,
                         mapq = bam[[1]]$mapq,
                         cigar = bam[[1]]$cigar,
                         mpos = bam[[1]]$mpos,
                         isize = bam[[1]]$isize)
  
  reads_dt[, bam_file := bam_path]
  return(reads_dt)
}

# get alleles
fasta <- read.fasta(reference_fasta_path)
hla_alleles <- names(fasta)
patient_genes <- unique(gsub("(hla_.*?)_.*", "\\1", hla_alleles))

reads <- data.table()
for(bam_path in hla_bam_paths){
  reads <- rbindlist(list(reads, 
                          make_read_table(bam_path)), use.name = TRUE)
}
reads[, gene := gsub("(hla_.*?)_.*", "\\1", chrom)]

# double check that a read and it's mate always maps to the same allele
if(nrow(reads[chrom != mate_chrom]) > 0){
  stop("A read and it's mate map to different alleles. Something has gone wrong as we filtered for only properly paired reads when making the HLA BAM file.")
}

# get reads that map to more than one gene
gene_dt <- unique(reads[,c("read", "gene")])
reads_mapping_multi_gene_dt <- gene_dt[,.N,by = c("read")][N>1]

gene_multimapping_reads_dt <- gene_dt[read %in% reads_mapping_multi_gene_dt$read]
n_gene_multi_mapping_reads <- gene_multimapping_reads_dt[, .N,by = "gene"]
setnames(n_gene_multi_mapping_reads, "N", "n_multi_gene_reads")

# remove reads that map to more than one gene
reads_single_gene <- reads[!read %in% reads_mapping_multi_gene_dt$read]
allele_read_count <- data.table()
for(g in patient_genes){
  
  gene_alleles <- grep(g, hla_alleles, value = TRUE)
  
  x <- unique(reads_single_gene[gene == g, c("read", "chrom")])
  read_count <- x[,.N,by = c("read")]
  reads_mapping_two_alleles <- unique(read_count[N==2]$read)
  
  reads_mapping_single_allele <- unique(x[!read %in% reads_mapping_two_alleles])
  
  # double check
  if(nrow(reads_mapping_single_allele[,.N,by = "read"][N!=1]) > 0){
    stop("This should only have reads mapping to single alleles")
  }
  
  n_bams <- length(unique(reads_mapping_single_allele$chrom))
  n_alleles <- length(gene_alleles)
  
  if(!n_alleles %in% c(1,2)){
    stop("There should be one or two alleles per gene")
  }
  
  allele_count <- reads_mapping_single_allele[,.N,by = "chrom"]
  setnames(allele_count, c("allele", "n_reads_map_uniquely"))
  allele_count[, sample_name := sample_name]
  
  if(length(gene_alleles) == 2){
    # heterozygous gene
    if(n_bams == 2){
      
      # two alleles with BAM
      allele_count[, allele_num := c("allele1", "allele2")]
      y <- dcast(allele_count, sample_name ~ allele_num, value.var = "n_reads_map_uniquely")
      setnames(y, 
               c("allele1", "allele2"),
               c("allele1_n_reads_mapping_uniquely",
                 "allele2_n_reads_mapping_uniquely"))
      y[,allele1 := allele_count[allele_num == "allele1"]$allele]
      y[,allele2 := allele_count[allele_num == "allele2"]$allele]
      y[,n_multi_allele_reads := length(reads_mapping_two_alleles)]  
      y[, homozygous := FALSE]
    }else if(n_bams == 1){
      # only one allele has filtered bam
      setnames(allele_count, c("allele1", "allele1_n_reads_mapping_uniquely",
                               "sample_name"))
      
      missing_allele <- gene_alleles[!gene_alleles %in% allele_count$allele1]
      allele_count[, allele2 := missing_allele ]
      allele_count[, allele2_n_reads_mapping_uniquely := 0]
      allele_count[,n_multi_allele_reads := 0]
      y <- allele_count
      y[, homozygous := FALSE]
    }else if(n_bams == 0){
      y <- data.table(sample_name = sample_name,
                      allele1 = gene_alleles[1],
                      allele2 = gene_alleles[2],
                      allele1_n_reads_mapping_uniquely = 0,
                      allele2_n_reads_mapping_uniquely = 0,
                      n_multi_allele_reads = 0,
                      homozygous = FALSE)
    }else{
      stop("With a heterozygous gene you should only have 0, 1 or 2 BAM files")
    }
    
  }else{
    # homozygous
    if(n_bams == 1){
      allele_count[, homozygous := TRUE]
      setnames(allele_count, c("allele1", "allele1_n_reads_mapping_uniquely",
                               "sample_name", "homozygous"))
      allele_count[, allele2 := as.character(NA)]
      allele_count[, allele2_n_reads_mapping_uniquely := as.numeric(NA)]
      allele_count[,n_multi_allele_reads := 0]
      y <- allele_count
    }else if(n_bams == 0){
      y <- data.table(sample_name = sample_name,
                      allele1 = gene_alleles,
                      allele2 = as.character(NA),
                      allele1_n_reads_mapping_uniquely = 0,
                      allele2_n_reads_mapping_uniquely = as.numeric(NA),
                      n_multi_allele_reads = 0,
                      homozygous = TRUE)
    }else{
      stop("With a homozygous gene you should only have 0 or 1 BAM files")
    }
  }
  
  allele_read_count <- rbindlist(list(allele_read_count, y), use.names = TRUE)
}

# add in the number of reads that map to multigene
allele_read_count[, gene := toupper(gsub(".*hla_(.*?)_.*", "\\1", allele1))]
n_gene_multi_mapping_reads[, gene := toupper(gsub("hla_(.*)", "\\1", gene))]
allele_read_count <- merge(allele_read_count, n_gene_multi_mapping_reads, by = "gene", all.x = TRUE)
allele_read_count[is.na(n_multi_gene_reads), n_multi_gene_reads := 0]
allele_read_count[, patient := patient]

fwrite(unique(reads_single_gene[,"read"]), file = paste0(out_file_prefix, "_reads_mapping_one_allele.csv"))
fwrite(allele_read_count, file = paste0(out_file_prefix, "_allele_read_count.csv"))
fwrite(gene_multimapping_reads_dt, file = paste0(out_file_prefix, "_gene_multimapping_reads.csv"))
