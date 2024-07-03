check_sj_within_feature <- function(sj_chr, sj_start, sj_end, gtf){
  
  sj_in_feature_dt <- gtf[seqname == sj_chr & 
                            sj_start > feature_start &
                            sj_end < feature_end]
  if(nrow(sj_in_feature_dt) > 0){
    
    if(nrow(sj_in_feature_dt) != 1){
      stop("what?")
    }
    
    output <- list(sj_in_feature = TRUE,
                   sj_in_feature_name = sj_in_feature_dt$feature_name)
    
    
  }else{
    output <- list(sj_in_feature = FALSE,
                   sj_in_feature_name = as.character(NA))
  }
  return(output)
}

check_exon_skip <- function(sj_chr, sj_start, sj_end, gtf){
  
  full_exon_skip_dt <- gtf[type == "exon" & seqname == sj_chr & 
                             sj_start <= feature_start &
                             sj_end >= feature_end]
  
  if(nrow(full_exon_skip_dt) > 0){
    
    output <- list(full_exon_skip = TRUE,
                   n_exon_skipped = nrow(full_exon_skip_dt),
                   exon_skipped_names = paste0(full_exon_skip_dt$feature_name, collapse = ";"))
    
  }else{
    output <- list(full_exon_skip = FALSE,
                   n_exon_skipped = as.numeric(NA),
                   exon_skipped_names = as.character(NA))
  }
  return(output)
}

check_exon_end_skip_forward_strand <- function(sj_chr, sj_start, gtf){
  
  exon_end_skipped_dt <- gtf[type == "exon" & seqname == sj_chr & 
                               feature_start < sj_start &
                               sj_start < feature_end]
  
  if(nrow(exon_end_skipped_dt) > 0){
    
    if(nrow(exon_end_skipped_dt) != 1){
      stop("Can only have one end skipped?")
    }
    
    # check start end and how many bases skipped    
    # get the exon start and end
    exon_start <- gtf[seqname == sj_chr & feature_name == exon_end_skipped_dt$feature_name]$feature_start
    exon_end <- gtf[seqname == sj_chr & feature_name == exon_end_skipped_dt$feature_name]$feature_end
    
    output <- list(end_exon_skip = TRUE,
                   end_exon_skipped_name = exon_end_skipped_dt$feature_name,
                   end_exon_skipped_start = exon_start,
                   end_exon_skipped_end = exon_end)
    
  }else{
    
    output <- list(end_exon_skip = FALSE,
                   end_exon_skipped_name = as.character(NA),
                   end_exon_skipped_start = as.numeric(NA),
                   end_exon_skipped_end = as.numeric(NA))
  }
  return(output)
}

check_exon_start_skip_forward_strand <- function(sj_chr, sj_end, gtf){
  
  exon_start_skipped_dt <- gtf[type == "exon" & seqname == sj_chr & 
                                 feature_start < sj_end &
                                 sj_end < feature_end]
  
  if(nrow(exon_start_skipped_dt) > 0){
    
    if(nrow(exon_start_skipped_dt) != 1){
      stop("Can only have one end skipped?")
    }
    
    # check start end and how many bases skipped    
    # get the exon start and end
    exon_start <- gtf[seqname == sj_chr & feature_name == exon_start_skipped_dt$feature_name]$feature_start
    exon_end <- gtf[seqname == sj_chr & feature_name == exon_start_skipped_dt$feature_name]$feature_end
    
    output <- list(start_exon_skip = TRUE,
                   start_exon_skipped_name = exon_start_skipped_dt$feature_name,
                   start_exon_skipped_start = exon_start,
                   start_exon_skipped_end = exon_end)
    
  }else{
    
    output <- list(start_exon_skip = FALSE,
                   start_exon_skipped_name = as.character(NA),
                   start_exon_skipped_start = as.numeric(NA),
                   start_exon_skipped_end = as.numeric(NA))
  }
  return(output)
}

check_exon_end_skip_reverse_strand <- function(sj_chr, sj_end, gtf){
  
  # if the gene is on the reverse strand, then "feature_end" will be the start
  # and "feature_start" will be the end
  exon_end_skipped_dt <- gtf[type == "exon" & seqname == sj_chr & 
                               feature_start < sj_end &
                               sj_end < feature_end]
  
  if(nrow(exon_end_skipped_dt) > 0){
    
    if(nrow(exon_end_skipped_dt) != 1){
      stop("Can only have one end skipped?")
    }
    
    # check start end and how many bases skipped    
    # get the exon start and end
    exon_start <- gtf[seqname == sj_chr & feature_name == exon_end_skipped_dt$feature_name]$feature_end
    exon_end <- gtf[seqname == sj_chr & feature_name == exon_end_skipped_dt$feature_name]$feature_start
    
    output <- list(end_exon_skip = TRUE,
                   end_exon_skipped_name = exon_end_skipped_dt$feature_name,
                   end_exon_skipped_start = exon_start,
                   end_exon_skipped_end = exon_end)
    
  }else{
    
    output <- list(end_exon_skip = FALSE,
                   end_exon_skipped_name = as.character(NA),
                   end_exon_skipped_start = as.numeric(NA),
                   end_exon_skipped_end = as.numeric(NA))
  }
  return(output)
}

check_exon_start_skip_reverse_strand <- function(sj_chr, sj_start, gtf){
  
  # if the gene is on the reverse strand, then "feature_end" will be the start
  # and "feature_start" will be the end
  exon_start_skipped_dt <- gtf[type == "exon" & seqname == sj_chr & 
                                 feature_start < sj_start &
                                 sj_start < feature_end]
  
  if(nrow(exon_start_skipped_dt) > 0){
    
    if(nrow(exon_start_skipped_dt) != 1){
      stop("Can only have one end skipped?")
    }
    
    # check start end and how many bases skipped    
    # get the exon start and end
    exon_start <- gtf[seqname == sj_chr & feature_name == exon_start_skipped_dt$feature_name]$feature_end
    exon_end <- gtf[seqname == sj_chr & feature_name == exon_start_skipped_dt$feature_name]$feature_start
    
    output <- list(start_exon_skip = TRUE,
                   start_exon_skipped_name = exon_start_skipped_dt$feature_name,
                   start_exon_skipped_start = exon_start,
                   start_exon_skipped_end = exon_end)
    
  }else{
    
    output <- list(start_exon_skip = FALSE,
                   start_exon_skipped_name = as.character(NA),
                   start_exon_skipped_start = as.numeric(NA),
                   start_exon_skipped_end = as.numeric(NA))
  }
  return(output)
}

check_intron_start_retained_forward_strand <- function(sj_chr, sj_start, gtf){
  
  intron_start_retained_dt <- gtf[type == "intron" & seqname == sj_chr & 
                                    feature_start < sj_start &
                                    sj_start < feature_end]
  
  if(nrow(intron_start_retained_dt) > 0){
    
    if(nrow(intron_start_retained_dt) != 1){
      stop("Can only have one start retained?")
    }
    
    # get the exon start and end
    intron_start <- gtf[seqname == sj_chr & feature_name == intron_start_retained_dt$feature_name]$feature_start
    intron_end <- gtf[seqname == sj_chr & feature_name == intron_start_retained_dt$feature_name]$feature_end
    
    output <- list(start_intron_retained = TRUE,
                   start_intron_retained_name = intron_start_retained_dt$feature_name,
                   start_intron_retained_start = intron_start,
                   start_intron_retained_end = intron_end)
    
  }else{
    
    output <- list(start_intron_retained = FALSE,
                   start_intron_retained_name = as.character(NA),
                   start_intron_retained_start = as.numeric(NA),
                   start_intron_retained_end = as.numeric(NA))
    
  }
  return(output)
}

check_intron_end_retained_forward_strand <- function(sj_chr, sj_end, gtf){
  
  intron_end_retained_dt <- gtf[type == "intron" & seqname == sj_chr & 
                                  feature_start < sj_end &
                                  sj_end < feature_end]
  
  if(nrow(intron_end_retained_dt) > 0){
    
    if(nrow(intron_end_retained_dt) != 1){
      stop("Can only have one start retained?")
    }
    
    # get the exon start and end
    intron_start <- gtf[seqname == sj_chr & feature_name == intron_end_retained_dt$feature_name]$feature_start
    intron_end <- gtf[seqname == sj_chr & feature_name == intron_end_retained_dt$feature_name]$feature_end
    
    output <- list(end_intron_retained = TRUE,
                   end_intron_retained_name = intron_end_retained_dt$feature_name,
                   end_intron_retained_start = intron_start,
                   end_intron_retained_end = intron_end)
    
  }else{
    
    output <- list(end_intron_retained = FALSE,
                   end_intron_retained_name = as.character(NA),
                   end_intron_retained_start = as.numeric(NA),
                   end_intron_retained_end = as.numeric(NA))
    
  }
  return(output)
}

check_intron_start_retained_reverse_strand <- function(sj_chr, sj_end, gtf){
  
  intron_start_retained_dt <- gtf[type == "intron" & seqname == sj_chr & 
                                    feature_start < sj_end &
                                    sj_end < feature_end]
  
  if(nrow(intron_start_retained_dt) > 0){
    
    if(nrow(intron_start_retained_dt) != 1){
      stop("Can only have one start retained?")
    }
    
    # get the exon start and end
    intron_start <- gtf[seqname == sj_chr & feature_name == intron_start_retained_dt$feature_name]$feature_end
    intron_end <- gtf[seqname == sj_chr & feature_name == intron_start_retained_dt$feature_name]$feature_start
    
    output <- list(start_intron_retained = TRUE,
                   start_intron_retained_name = intron_start_retained_dt$feature_name,
                   start_intron_retained_start = intron_start,
                   start_intron_retained_end = intron_end)
    
  }else{
    
    output <- list(start_intron_retained = FALSE,
                   start_intron_retained_name = as.character(NA),
                   start_intron_retained_start = as.numeric(NA),
                   start_intron_retained_end = as.numeric(NA))
    
  }
  return(output)
}

check_intron_end_retained_reverse_strand <- function(sj_chr, sj_start, gtf){
  
  intron_end_retained_dt <- gtf[type == "intron" & seqname == sj_chr & 
                                  feature_start < sj_start &
                                  sj_start < feature_end]
  
  if(nrow(intron_end_retained_dt) > 0){
    
    if(nrow(intron_end_retained_dt) != 1){
      stop("Can only have one start retained?")
    }
    
    # get the exon start and end
    intron_start <- gtf[seqname == sj_chr & feature_name == intron_end_retained_dt$feature_name]$feature_end
    intron_end <- gtf[seqname == sj_chr & feature_name == intron_end_retained_dt$feature_name]$feature_start
    
    output <- list(end_intron_retained = TRUE,
                   end_intron_retained_name = intron_end_retained_dt$feature_name,
                   end_intron_retained_start = intron_start,
                   end_intron_retained_end = intron_end)
    
  }else{
    
    output <- list(end_intron_retained = FALSE,
                   end_intron_retained_name = as.character(NA),
                   end_intron_retained_start = as.numeric(NA),
                   end_intron_retained_end = as.numeric(NA))
    
  }
  return(output)
}

get_codon_tab <- function(nuc_seq, codon_mapping){
  
  codons <- c()
  idx <- 1
  for(i in 1:(length(nuc_seq)/3)){
    codon <- paste0(nuc_seq[idx:(idx+2)], collapse = "")
    idx <- idx + 3
    codons <- c(codons, codon)
  }
  
  codon_dt <- data.table(codon = codons)
  codon_dt[, pos := 1:nrow(codon_dt)]
  codon_dt <- merge(codon_dt, codon_aa_tab, by = "codon", all.x = TRUE)
  codon_dt <- codon_dt[order(pos)]
  
  return(codon_dt)
}

