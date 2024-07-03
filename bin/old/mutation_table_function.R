
get_bam_read_count <- function(bam_path, chr, start, stop, ref, alt, mutation_type){
  
  if(mutation_type == "SNV"){
    
    wxs_reads_output <- as.data.table(bam2R(file = bam_path,
                                            chr = chr,
                                            start = start,
                                            stop = start))
    
    bam_ref_count <- as.numeric(get_snp_allele_count(ref, wxs_reads_output))
    bam_alt_count <- as.numeric(get_snp_allele_count(alt, wxs_reads_output))
    bam_N_count <- as.numeric(get_snp_allele_count("N", wxs_reads_output))
    
  }
  
  if(mutation_type == "MNV"){
    
    wxs_reads_output <- as.data.table(bam2R(file = bam_path,
                                            chr = chr,
                                            start = start,
                                            stop = stop))
    
    bam_ref_count <- as.numeric(get_mnv_allele_count(ref, wxs_reads_output)$allele_count)
    bam_alt_count <- as.numeric(get_mnv_allele_count(alt, wxs_reads_output)$allele_count)
    bam_N_count <- as.numeric(get_mnv_allele_count(ref, wxs_reads_output)$n_count)
    
  }
  
  if(mutation_type == c("INS")){
    
    wxs_reads_output <- as.data.table(bam2R(file = bam_path,
                                            chr = chr,
                                            start = start,
                                            stop = start))
    
    bam_ref_count <- as.numeric(get_ins_allele_count(ref, wxs_reads_output)$ref_count)
    bam_alt_count <- as.numeric(get_ins_allele_count(ref, wxs_reads_output)$alt_count)
    bam_N_count <- as.numeric(get_ins_allele_count(ref, wxs_reads_output)$n_count)
    
  }
  
  if(mutation_type == c("DEL")){
    
    wxs_reads_output <- as.data.table(bam2R(file = bam_path,
                                            chr = chr,
                                            start = start + 1,
                                            stop = stop))
    
    bam_ref_count <- as.numeric(get_del_allele_count(ref, wxs_reads_output)$ref_count)
    bam_alt_count <- as.numeric(get_del_allele_count(ref, wxs_reads_output)$alt_count)
    bam_N_count <- as.numeric(get_del_allele_count(ref, wxs_reads_output)$n_count)
    
  }
  
  return(list(bam_ref_count = bam_ref_count,
              bam_alt_count = bam_alt_count,
              bam_N_count = bam_N_count))
  
}

get_snp_allele_count <- function(allele, bam2R_dt){
  
  if(allele == "A"){
    allele_count <- bam2R_dt$A + bam2R_dt$a
  }
  
  if(allele == "C"){
    allele_count <- bam2R_dt$C + bam2R_dt$c
  }
  
  if(allele == "G"){
    allele_count <- bam2R_dt$G + bam2R_dt$g
  }
  
  if(allele == "T"){
    allele_count <- bam2R_dt$`T` + bam2R_dt$t
  }
  
  if(allele == "N"){
    allele_count <- bam2R_dt$`N` + bam2R_dt$n
  }
  
  return(allele_count)
}

get_mnv_allele_count <- function(allele, bam2R_dt){
  
  if(nchar(allele) != nrow(bam2R_dt)){
    stop("?")
  }
  
  split_allele <- strsplit(allele, split = "")[[1]]
  allele_counts <- c()
  row_idx <- 1
  for(a in split_allele){
    if(a == "A"){
      allele_counts <- c(allele_counts, bam2R_dt$A[row_idx] + bam2R_dt$a[row_idx])
    }
    
    if(a == "C"){
      allele_counts <- c(allele_counts, bam2R_dt$C[row_idx] + bam2R_dt$c[row_idx])
    }
    
    if(a == "G"){
      allele_counts <- c(allele_counts, bam2R_dt$G[row_idx] + bam2R_dt$g[row_idx])
    }
    
    if(a == "T"){
      allele_counts <- c(allele_counts, bam2R_dt$`T`[row_idx] + bam2R_dt$t[row_idx])
    }
    row_idx <- row_idx + 1
  }
  return(list(allele_count = mean(allele_counts),
              n_count = mean(bam2R_dt$N) + mean(bam2R_dt$n)))
}

get_del_allele_count <- function(ref, bam2R_dt){
  
  if((nchar(ref) - 1) != nrow(bam2R_dt)){
    stop("?")
  }
  
  alt_count <- mean(bam2R_dt$DEL) + mean(bam2R_dt$del)
  
  split_allele <- strsplit(ref, split = "")[[1]]
  allele_counts <- c()
  row_idx <- 1
  for(a in split_allele[2:length(split_allele)]){
    if(a == "A"){
      allele_counts <- c(allele_counts, bam2R_dt$A[row_idx] + bam2R_dt$a[row_idx])
    }
    
    if(a == "C"){
      allele_counts <- c(allele_counts, bam2R_dt$C[row_idx] + bam2R_dt$c[row_idx])
    }
    
    if(a == "G"){
      allele_counts <- c(allele_counts, bam2R_dt$G[row_idx] + bam2R_dt$g[row_idx])
    }
    
    if(a == "T"){
      allele_counts <- c(allele_counts, bam2R_dt$`T`[row_idx] + bam2R_dt$t[row_idx])
    }
    row_idx <- row_idx + 1
  }
  return(list(ref_count = mean(allele_counts),
              n_count = mean(bam2R_dt$N) + mean(bam2R_dt$n),
              alt_count = alt_count))
}

get_ins_allele_count <- function(ref, bam2R_dt){
  
  if(nrow(bam2R_dt) != 1){
    stop("The input bam2R_dt should be one row")
  }
  
  if(nchar(ref) != 1){
    stop("With insertions, REF should be length 1")
  }
  
  alt_count <- bam2R_dt$INS + bam2R_dt$ins
  
  if(ref == "A"){
    ref_count <- bam2R_dt$A + bam2R_dt$a
  }
  
  if(ref == "C"){
    ref_count <- bam2R_dt$C + bam2R_dt$c
  }
  
  if(ref == "G"){
    ref_count <- bam2R_dt$G + bam2R_dt$g
  }
  
  if(ref == "T"){
    ref_count <- bam2R_dt$`T` + bam2R_dt$t
  }
  
  return(list(ref_count = ref_count,
              alt_count =  alt_count,
              n_count = bam2R_dt$N + bam2R_dt$n))
}