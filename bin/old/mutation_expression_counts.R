suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(deepSNV))

args <- commandArgs(trailingOnly = TRUE)
rna_sample_name_hash <- args[1]
dna_sample_name_hash <- args[2]

mut_table_path <- list.files(pattern = "mut_table.csv")

x <- read.csv(mut_table_path, stringsAsFactors=F)
x$rna_sample_name_hash <- rna_sample_name_hash

for(x_idx in 1:nrow(x)){

a_name <- x[x_idx, "chr"]

# add the rna bam path
rna_bam_path <-  paste0(rna_sample_name_hash, ".", a_name,".sorted.filtered.bam")

if(length(rna_bam_path) > 1){
    stop("should be one rna bam")
}else if(length(rna_bam_path) == 1 & !is.na(rna_bam_path) & file.exists(rna_bam_path)){
    x[x_idx, "rna_bam"] <- rna_bam_path
}
   
# get the number of reads from the rna bam
    if(x[x_idx,"snv"] & length(rna_bam_path) == 1 & !is.na(rna_bam_path) & !is.na(x[x_idx,"pos_in_cds"])){
        rna_reads_output <- as.data.frame(bam2R(file = x[x_idx,"rna_bam"],
                                                chr = x[x_idx,"chr"],
                                                start = x[x_idx,"pos_in_cds"],
                                                stop = x[x_idx,"pos_in_cds"]))
        rna_reads_output <- rna_reads_output[,c("A", "C", "G", "T", "a", "c", "g", "t")]
        
        x[x_idx, "rna_A"] <- rna_reads_output$A
        x[x_idx, "rna_C"] <- rna_reads_output$C
        x[x_idx, "rna_G"] <- rna_reads_output$G
        x[x_idx, "rna_T"] <- rna_reads_output$`T`
        x[x_idx, "rna_a"] <- rna_reads_output$a
        x[x_idx, "rna_c"] <- rna_reads_output$c
        x[x_idx, "rna_g"] <- rna_reads_output$g
        x[x_idx, "rna_t"] <- rna_reads_output$t
    }
}

# overwrite csv
write.csv(x, file = paste0(dna_sample_name_hash, ".mut_table.csv"))
