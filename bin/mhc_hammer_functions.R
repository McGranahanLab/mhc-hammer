suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbeeswarm))

# from https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

calculate_cn_with_baf <- function(logr,baf,ploidy,purity,gamma){
  
  nA <- (purity-1+baf*2^(logr/gamma)*((1-purity)*2+purity*ploidy))/purity
  nB <- (purity-1-(baf-1)*2^(logr/gamma)*((1-purity)*2+purity*ploidy))/purity
  return(list(nA = nA, nB = nB))
  
}

calculate_logr_aib <- function(gl_library_size, 
                               tumour_library_size,
                               allele1_name,
                               allele2_name,
                               allele1_snp_bed_file,
                               allele2_snp_bed_file,
                               allele1_gl_reads_count_once_coverage,
                               allele1_tumour_reads_count_once_coverage,
                               allele2_gl_reads_count_once_coverage,
                               allele2_tumour_reads_count_once_coverage,
                               make_plot){
  
  mult_factor <- gl_library_size/tumour_library_size
  
  allele1_snps <- fread(allele1_snp_bed_file, sep = "\t")
  allele2_snps <- fread(allele2_snp_bed_file, sep = "\t")
  
  if(nrow(allele1_snps) == 0 | nrow(allele2_snps) == 0){
    return(list(paired_t_test = as.numeric(NA),
                paired_wilcoxon_test = as.numeric(NA),
                n_snps = 0,
                boxplot = ggplot()))
  }
  
  setnames(allele1_snps, c("CHROM", "chromStart", "chromEnd"))
  setnames(allele2_snps, c("CHROM", "chromStart", "chromEnd"))
  
  allele1_gl_coverage <- fread(allele1_gl_reads_count_once_coverage)
  setnames(allele1_gl_coverage, c("position", "allele1_gl_depth"))
  allele1_tumour_coverage <-fread(allele1_tumour_reads_count_once_coverage)
  setnames(allele1_tumour_coverage, c("position", "allele1_tumour_depth"))
  
  allele2_gl_coverage <- fread(allele2_gl_reads_count_once_coverage)
  setnames(allele2_gl_coverage, c("position", "allele2_gl_depth"))
  allele2_tumour_coverage <-fread(allele2_tumour_reads_count_once_coverage)
  setnames(allele2_tumour_coverage, c("position", "allele2_tumour_depth"))
  
  if(nrow(allele1_gl_coverage) == 0 | nrow(allele2_gl_coverage) == 0 |
     nrow(allele1_tumour_coverage) == 0 | nrow(allele2_tumour_coverage) == 0){
    # return(list(paired_t_test = as.numeric(NA),
    #             paired_wilcoxon_test = as.numeric(NA),
    #             n_snps = 0,
    #             boxplot = ggplot()))
    stop("Why are there no lines?")
  }
  
  coverage_table_at_mismatch <- data.table(allele1_mismatch_pos = allele1_snps$chromEnd,
                                           allele2_mismatch_pos = allele2_snps$chromEnd)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele1_gl_coverage,
                                      by.x = "allele1_mismatch_pos", by.y = "position", 
                                      all.x = TRUE)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele1_tumour_coverage,
                                      by.x = "allele1_mismatch_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele2_gl_coverage,
                                      by.x = "allele2_mismatch_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele2_tumour_coverage,
                                      by.x = "allele2_mismatch_pos", by.y = "position", all.x = TRUE)
  
  if(nrow(coverage_table_at_mismatch[is.na(allele1_gl_depth) | is.na(allele1_tumour_depth) |
                                     is.na(allele2_gl_depth) | is.na(allele2_tumour_depth)]) > 0){
    stop("why are there missing lines?")
  }
  
  coverage_table_at_mismatch[,allele1_tumour_depth := as.numeric(allele1_tumour_depth)]
  coverage_table_at_mismatch[,allele2_tumour_depth := as.numeric(allele2_tumour_depth)]
  
  coverage_table_at_mismatch[,allele1_gl_depth := as.numeric(allele1_gl_depth)]
  coverage_table_at_mismatch[,allele2_gl_depth := as.numeric(allele2_gl_depth)]
  
  coverage_table_at_mismatch[, allele1_logr := log2( (allele1_tumour_depth/allele1_gl_depth + 0.0001 ) * mult_factor )]
  coverage_table_at_mismatch[, allele2_logr := log2( (allele2_tumour_depth/allele2_gl_depth + 0.0001 ) * mult_factor )]
  
  # if one of the alleles has duplicated mismatches, the logr for the other allele
  # will change, so we take the median
  coverage_table_at_mismatch[, allele1_logr_no_dup := median(allele1_logr), by = allele2_mismatch_pos]
  coverage_table_at_mismatch[, allele2_logr_no_dup := median(allele2_logr), by = allele1_mismatch_pos]
  
  # mark and remove duplicated values
  coverage_table_at_mismatch[, allele1_dup_id := 1:.N, allele1_mismatch_pos]
  coverage_table_at_mismatch[, allele2_dup_id := 1:.N, allele2_mismatch_pos]
  coverage_table_at_mismatch <- coverage_table_at_mismatch[allele1_dup_id == 1 & allele2_dup_id == 1]
  
  n_snps <- nrow(coverage_table_at_mismatch)
  
  # get p values
  paired_t_test <- myTryCatch(t.test(coverage_table_at_mismatch[,allele1_logr_no_dup],
                                     coverage_table_at_mismatch[,allele2_logr_no_dup],
                                     paired=TRUE))
  
  paired_wilcoxon_test <- myTryCatch(wilcox.test(x = coverage_table_at_mismatch[,allele1_logr_no_dup],
                                                 y = coverage_table_at_mismatch[,allele2_logr_no_dup],
                                                 paired=TRUE))
  
  if(is.null(paired_t_test$value)){
    dna_aib_t_test <- NA
  }else{
    dna_aib_t_test <- paired_t_test$value$p.value
  }
  
  if(is.null(paired_wilcoxon_test$value)){
    dna_aib_wilcox_test <- NA
  }else{
    dna_aib_wilcox_test <- paired_wilcoxon_test$value$p.value
  }
  
  if(make_plot){
    
    # plot the logR
    allele1_logr_plot <- coverage_table_at_mismatch[,c("allele1_mismatch_pos",
                                                       "allele1_logr_no_dup")]
    allele1_logr_plot[,allele_name := allele1_name]
    allele1_logr_plot[,allele_col := "#377EB8"]
    
    allele2_logr_plot <- coverage_table_at_mismatch[,c("allele1_mismatch_pos",
                                                       "allele2_logr_no_dup")]
    allele2_logr_plot[,allele_name := allele2_name]
    allele2_logr_plot[,allele_col := "#E41A1C"]
    
    logr_plot_dat <- rbindlist(list(allele1_logr_plot, 
                                    allele2_logr_plot),
                               use.names = FALSE)
    
    col_dt <- unique(logr_plot_dat[,c("allele_name", "allele_col")])
    
    logr_boxplot_plot <- ggplot(logr_plot_dat, 
                                aes(allele_name, 
                                    allele1_logr_no_dup,
                                    colour = allele_name)) +
      geom_boxplot() +
      geom_quasirandom() +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.position = "none") +
      ylab("LogR") +
      ggtitle("DNA: LogR", 
              subtitle = paste0("t test p value=", 
                                formatC(paired_t_test$value$p.value, format = "e", digits = 2),
                                "\nwilcoxon p value=", 
                                formatC(paired_wilcoxon_test$value$p.value, format = "e", digits = 2))) +
      scale_colour_manual(breaks = c(col_dt$allele_name),
                          values = c(col_dt$allele_col))
  }else{
    logr_boxplot_plot <-  ggplot()
  }
  
  return(list(paired_t_test = dna_aib_t_test,
              paired_wilcoxon_test = dna_aib_wilcox_test,
              n_snps = n_snps,
              boxplot = logr_boxplot_plot))
  
}

calculate_cn <- function(bin_size =  150, gamma = 1,
                         allele1_name,
                         allele2_name,
                         allele1_gl_coverage_file,
                         allele1_tumour_coverage_file,
                         allele2_gl_coverage_file,
                         allele2_tumour_coverage_file,
                         allele1_snp_bed_file,
                         allele2_snp_bed_file,
                         purity,
                         ploidy,
                         gl_library_size,
                         tumour_library_size,
                         gtf_path,
                         make_plot){
  
  allele1_snp <- fread(allele1_snp_bed_file, sep = "\t")
  allele2_snp <- fread(allele2_snp_bed_file, sep = "\t")
  
  if(nrow(allele1_snp) == 0 | nrow(allele2_snp) == 0){
    return(list(cn1 = as.numeric(NA),
                cn1_lower = as.numeric(NA),
                cn1_upper = as.numeric(NA),
                cn2 = as.numeric(NA),
                cn2_lower = as.numeric(NA),
                cn2_upper = as.numeric(NA),
                cn1_binned = as.numeric(NA),
                cn1_binned_lower = as.numeric(NA),
                cn1_binned_upper = as.numeric(NA),
                cn2_binned  = as.numeric(NA),
                cn2_binned_lower = as.numeric(NA),
                cn2_binned_upper = as.numeric(NA),
                cn_boxplot_plot = ggplot(),
                cn_bin_boxplot_plot =  ggplot(),
                allele1_cn_plot = ggplot(),
                allele2_cn_plot = ggplot(),
                allele1_logr_plot = ggplot(),
                allele2_logr_plot = ggplot(),
                allele1_coverage_plot = ggplot(),
                allele2_coverage_plot = ggplot(),
                n_snps = 0,
                n_bins = as.numeric(NA)))
  }
  
  setnames(allele1_snp, c("chr", "start", "end"))
  setnames(allele2_snp, c("chr", "start", "end"))
  
  allele1_gl_coverage <- fread(allele1_gl_coverage_file)
  setnames(allele1_gl_coverage, c("position", "allele1_gl_depth"))
  allele1_tumour_coverage <-fread(allele1_tumour_coverage_file)
  setnames(allele1_tumour_coverage, c("position", "allele1_tumour_depth"))
  
  allele2_gl_coverage <- fread(allele2_gl_coverage_file)
  setnames(allele2_gl_coverage, c("position", "allele2_gl_depth"))
  allele2_tumour_coverage <-fread(allele2_tumour_coverage_file)
  setnames(allele2_tumour_coverage, c("position", "allele2_tumour_depth"))
  
  mult_factor <- gl_library_size/tumour_library_size
  
  # create a table with the coverage at the mismatch positions
  coverage_table_at_snps <- data.table(allele1_snp_pos = allele1_snp$end,
                                       allele2_snp_pos = allele2_snp$end)
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele1_gl_coverage,
                                  by.x = "allele1_snp_pos", by.y = "position", 
                                  all.x = TRUE)
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele1_tumour_coverage,
                                  by.x = "allele1_snp_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele2_gl_coverage,
                                  by.x = "allele2_snp_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele2_tumour_coverage,
                                  by.x = "allele2_snp_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_snps[, allele1_log2 := log2( allele1_tumour_depth/allele1_gl_depth * mult_factor )]
  coverage_table_at_snps[, allele2_log2 := log2( allele2_tumour_depth/allele2_gl_depth * mult_factor )]
  
  # if a snp position is repeated, eg if the other allele has a gap, then
  # the depth is repeated at the positions across the gap
  # we take the mean of the depth at the gaps is taken
  coverage_table_at_snps[, allele2_logr_no_dup := mean(allele2_log2), by = allele1_snp_pos]
  coverage_table_at_snps[, allele2_tumour_depth_no_dup := mean(allele2_tumour_depth), by = allele1_snp_pos]
  coverage_table_at_snps[, allele2_gl_depth_no_dup := mean(allele2_gl_depth), by = allele1_snp_pos]
  
  coverage_table_at_snps[, allele1_logr_no_dup := mean(allele1_log2), by = allele2_snp_pos]
  coverage_table_at_snps[, allele1_tumour_depth_no_dup := mean(allele1_tumour_depth), by = allele2_snp_pos]
  coverage_table_at_snps[, allele1_gl_depth_no_dup := mean(allele1_gl_depth), by = allele2_snp_pos]
  
  coverage_table_at_snps[,allele1_dup_id := 1:.N, allele1_snp_pos]
  coverage_table_at_snps[,allele2_dup_id := 1:.N, allele2_snp_pos]
  coverage_table_at_snps <- coverage_table_at_snps[allele1_dup_id ==1 & allele2_dup_id == 1]
  coverage_table_at_snps[,allele1_dup_id := NULL]
  coverage_table_at_snps[,allele2_dup_id := NULL]
  
  n_snps <- nrow(coverage_table_at_snps)
  
  coverage_table_at_snps[, baf := ( allele1_tumour_depth_no_dup /
                                      ( allele1_tumour_depth_no_dup + allele2_tumour_depth_no_dup ) )]
  
  coverage_table_at_snps[, logr := log2((( allele1_tumour_depth_no_dup + allele2_tumour_depth_no_dup ) /
                                           ( allele1_gl_depth_no_dup + allele2_gl_depth_no_dup )) * mult_factor)]
  
  # create a table with the coverage at the all filtered positions
  coverage_table <- merge(allele1_gl_coverage, allele1_tumour_coverage,
                          by = c("position"), all = TRUE)
  coverage_table <- merge(coverage_table, allele2_gl_coverage,
                          by = c("position"), all = TRUE)
  coverage_table <- merge(coverage_table, allele2_tumour_coverage,
                          by = c("position"), all = TRUE)
  
  coverage_table[, allele1_logr := log2(allele1_tumour_depth / allele1_gl_depth * mult_factor)]
  coverage_table[, allele2_logr := log2(allele2_tumour_depth / allele2_gl_depth * mult_factor)]
  
  # get binned values of the coverage table
  min_position <- min(c(allele1_tumour_coverage[,position], allele2_tumour_coverage[,position]))
  max_position   <- max(c(allele1_tumour_coverage[,position], allele2_tumour_coverage[,position]))
  bins <- seq(min_position, max_position, by = bin_size)
  bins <- c(bins[-length(bins)],max_position+1)
  
  n_bins <- length(bins)
  
  if(length(bins) <= 5){
    
    median_cn1_binned <- NA
    cn1_binned_lower <- NA
    cn1_binned_upper <- NA
    
    median_cn2_binned <- NA
    cn2_binned_lower <- NA
    cn2_binned_upper <- NA
    
    cn_bin_boxplot_plot <- ggplot()
    
  }else{
    coverage_table_binned <- data.table()
    for(i in 1:(length(bins) - 1)){
      
      lower_bound <- bins[i]
      upper_bound <- bins[i + 1]
      
      tmp_dt <- data.table(bin_start = lower_bound,
                           bin_end = upper_bound,
                           allele1_logr_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele1_logr, na.rm = TRUE ),
                           allele2_logr_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele2_logr, na.rm = TRUE ),
                           allele1_gl_median_depth_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele1_gl_depth, na.rm = TRUE ),
                           allele1_tumour_median_depth_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele1_tumour_depth, na.rm = TRUE ),
                           allele2_gl_median_depth_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele2_gl_depth, na.rm = TRUE ),
                           allele2_tumour_median_depth_bin = median( coverage_table[position %in% lower_bound:upper_bound]$allele2_tumour_depth, na.rm = TRUE ))
      
      coverage_table_binned <- rbindlist( list(coverage_table_binned, tmp_dt), use.names = TRUE )
      
      # add the bins to the depth table
      coverage_table_at_snps[lower_bound <= allele1_snp_pos & upper_bound > allele1_snp_pos,
                             c("bin_start", "bin_end") := .(lower_bound, upper_bound)]
    }
    
    coverage_table_binned[, tumour_depth_bin := allele1_tumour_median_depth_bin + allele2_tumour_median_depth_bin]
    coverage_table_binned[, gl_depth_bin := allele1_gl_median_depth_bin + allele2_gl_median_depth_bin]
    coverage_table_binned[, logr_bin := log2(tumour_depth_bin / gl_depth_bin * mult_factor)]
    
    coverage_table_at_snps <- merge(coverage_table_at_snps,
                                    coverage_table_binned,
                                    by = c("bin_start", "bin_end"),
                                    all.x = TRUE)
    # calculate the copy number using the binned log R
    allele_cn_binned <- calculate_cn_with_baf(logr = coverage_table_at_snps[,logr_bin],
                                              baf = coverage_table_at_snps[,baf],
                                              ploidy = ploidy,
                                              purity = purity,
                                              gamma = gamma)
    
    coverage_table_at_snps[, cn1_binned := allele_cn_binned$nA]
    coverage_table_at_snps[, cn2_binned := allele_cn_binned$nB]
    
    # only count each bin once
    coverage_table_at_snps[,bin_dup_id := 1:.N, bin_start]
    
    median_cn1_binned <- median(coverage_table_at_snps[bin_dup_id == 1,cn1_binned], na.rm = TRUE)
    median_cn2_binned <- median(coverage_table_at_snps[bin_dup_id == 1,cn2_binned], na.rm = TRUE)
    
    # get the confidence interval for cn1
    cn1_binned_t_test   <- myTryCatch(t.test(coverage_table_at_snps[bin_dup_id == 1, cn1_binned]))
    if(is.null(cn1_binned_t_test$value)){
      cn1_binned_lower <- NA
      cn1_binned_upper <- NA
    }else{
      cn1_binned_lower <- cn1_binned_t_test$value$conf.int[1]
      cn1_binned_upper <- cn1_binned_t_test$value$conf.int[2]
    }
    
    # get the confidence interval for cn2
    cn2_binned_t_test <- myTryCatch(t.test(coverage_table_at_snps[bin_dup_id == 1, cn2_binned]))
    if(is.null(cn2_binned_t_test$value)){
      cn2_binned_lower <- NA
      cn2_binned_upper <- NA
    }else{
      cn2_binned_lower <- cn2_binned_t_test$value$conf.int[1]
      cn2_binned_upper <- cn2_binned_t_test$value$conf.int[2]
    }
    
    if(make_plot){
      # plot the cn from the binned logr values
      allele1_cn_bin <- coverage_table_at_snps[bin_dup_id == 1,
                                               c("allele1_snp_pos",
                                                 "cn1_binned")]
      allele1_cn_bin[,allele_name := allele1_name]
      allele1_cn_bin[,allele_col := "#377EB8"]
      
      allele2_cn_bin <- coverage_table_at_snps[bin_dup_id == 1,
                                               c("allele1_snp_pos",
                                                 "cn2_binned")]
      
      allele2_cn_bin[,allele_name := allele2_name]
      allele2_cn_bin[,allele_col := "#E41A1C"]
      
      cn_bin_dat <- rbindlist(list(allele1_cn_bin,
                                   allele2_cn_bin),
                              use.names = FALSE)
      col_dt <- unique(cn_bin_dat[,c("allele_name", "allele_col")])
      
      cn_bin_boxplot_plot <- ggplot(cn_bin_dat,
                                    aes(allele_name,
                                        cn1_binned,
                                        colour = allele_name)) +
        geom_boxplot() +
        geom_quasirandom() +
        theme_bw() +
        theme(text = element_text(size= 20),
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
              legend.position = "none") +
        ylab("Copy number\n(using binned log R)") +
        ggtitle("DNA: Copy number (binned)") +
        scale_colour_manual(breaks = c(col_dt$allele_name),
                            values = c(col_dt$allele_col))
    }else{
      cn_bin_boxplot_plot <- ggplot()
    }
  }
  
  # calculate the copy number using the non binned log R
  allele_cn <- calculate_cn_with_baf(logr = coverage_table_at_snps[,logr],
                                     baf = coverage_table_at_snps[,baf],
                                     ploidy = ploidy,
                                     purity = purity,
                                     gamma = gamma)
  
  coverage_table_at_snps[, cn1 := allele_cn$nA]
  coverage_table_at_snps[, cn2 := allele_cn$nB]
  
  median_cn1 <- median(coverage_table_at_snps[,cn1], na.rm = TRUE)
  median_cn2 <- median(coverage_table_at_snps[,cn2], na.rm = TRUE)
  
  # get the confidence interval for cn1
  cn1_t_test   <- myTryCatch(t.test(coverage_table_at_snps[,cn1]))
  if(is.null(cn1_t_test$value)){
    cn1_lower <- NA
    cn1_upper <- NA
  }else{
    cn1_lower <- cn1_t_test$value$conf.int[1]
    cn1_upper <- cn1_t_test$value$conf.int[2]
  }
  
  # get the confidence interval for cn2
  cn2_t_test   <- myTryCatch(t.test(coverage_table_at_snps[,cn2]))
  if(is.null(cn2_t_test$value)){
    cn2_lower <- NA
    cn2_upper <- NA
  }else{
    cn2_lower <- cn2_t_test$value$conf.int[1]
    cn2_upper <- cn2_t_test$value$conf.int[2]
  }
  
  # plot the cn
  if(make_plot){
    allele1_cn <- coverage_table_at_snps[,c("allele1_snp_pos",
                                            "cn1")]
    allele1_cn[,allele_name := allele1_name]
    allele1_cn[,allele_col := "#377EB8"]
    
    allele2_cn <- coverage_table_at_snps[,c("allele1_snp_pos",
                                            "cn2")]
    
    allele2_cn[,allele_name := allele2_name]
    allele2_cn[,allele_col := "#E41A1C"]
    
    cn_dat <- rbindlist(list(allele1_cn,
                             allele2_cn),
                        use.names = FALSE)
    
    col_dt <- unique(cn_dat[,c("allele_name", "allele_col")])
    
    cn_boxplot_plot <- ggplot(cn_dat,
                              aes(allele_name,
                                  cn1,
                                  colour = allele_name)) +
      geom_boxplot() +
      geom_quasirandom() +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.position = "none") +
      ylab("Copy number") +
      ggtitle("DNA: Copy number") +
      scale_colour_manual(breaks = c(col_dt$allele_name),
                          values = c(col_dt$allele_col))
    
    
    gtf <- fread(gtf_path)
    
    allele1_tab_cn <- gtf[V1 == allele1_name & V3 == "exon"]
    allele1_tab_cn[,y_start := min(c(cn_dat[,cn1], 0), na.rm = TRUE)]
    allele1_tab_cn[,y_end := max(cn_dat[,cn1], na.rm = TRUE)]
    
    allele1_cn_plot <- ggplot() +
      geom_rect(data  = allele1_tab_cn,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_hline(yintercept = 0.5) +
      geom_hline(yintercept = median_cn1, colour  = col_dt[allele_name == allele1_name]$allele_col, linetype = "dashed") +
      geom_line(data = cn_dat[allele_name == allele1_name],
                aes(allele1_snp_pos,
                    cn1), colour  = col_dt[allele_name == allele1_name]$allele_col) +
      geom_point(data = cn_dat[allele_name == allele1_name],
                 aes(allele1_snp_pos,
                     cn1),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      xlab("Position in  allele") +
      ylab("CN") +
      ggtitle(allele1_name)
    
    allele2_tab_cn <- gtf[V1 == allele2_name & V3 == "exon"]
    allele2_tab_cn[,y_start := min(c(cn_dat$cn1, 0), na.rm = TRUE)]
    allele2_tab_cn[,y_end := max(cn_dat$cn1, na.rm = TRUE)]
    
    allele2_cn_plot <- ggplot() +
      geom_rect(data  = allele2_tab_cn,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_hline(yintercept = 0.5) +
      geom_hline(yintercept = median_cn2, colour  = col_dt[allele_name == allele2_name]$allele_col, linetype = "dashed") +
      geom_line(data = cn_dat[allele_name == allele2_name],
                aes(allele1_snp_pos,
                    cn1), colour  = col_dt[allele_name == allele2_name]$allele_col) +
      geom_point(data = cn_dat[allele_name == allele2_name],
                 aes(allele1_snp_pos,
                     cn1),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      xlab("Position in  allele") +
      ylab("CN") +
      ggtitle(allele2_name)
    
    # plot the logr
    allele1_logr <- coverage_table_at_snps[,c("allele1_snp_pos",
                                              "allele1_logr_no_dup")]
    allele1_logr[,allele_name := allele1_name]
    allele1_logr[,allele_col := "#377EB8"]
    
    allele2_logr <- coverage_table_at_snps[,c("allele1_snp_pos",
                                              "allele2_logr_no_dup")]
    
    allele2_logr[,allele_name := allele2_name]
    allele2_logr[,allele_col := "#E41A1C"]
    
    logr_dat <- rbindlist(list(allele1_logr,
                               allele2_logr),
                          use.names = FALSE)
    
    allele1_tab_logr <- gtf[V1 == allele1_name & V3 == "exon"]
    allele1_tab_logr[,y_start := min(logr_dat$allele1_logr_no_dup, na.rm = TRUE)]
    allele1_tab_logr[,y_end := max(logr_dat$allele1_logr_no_dup, na.rm = TRUE)]
    
    allele1_logr_plot <- ggplot() +
      geom_rect(data  = allele1_tab_logr,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_line(data = logr_dat[allele_name == allele1_name],
                aes(allele1_snp_pos,
                    allele1_logr_no_dup), colour  = col_dt[allele_name == allele1_name]$allele_col) +
      geom_point(data = logr_dat[allele_name == allele1_name],
                 aes(allele1_snp_pos,
                     allele1_logr_no_dup),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      xlab("Position in  allele") +
      ylab("logR") +
      ggtitle(allele1_name)
    
    allele2_tab_logr <- gtf[V1 == allele2_name & V3 == "exon"]
    allele2_tab_logr[,y_start := min(logr_dat$allele1_logr_no_dup, na.rm = TRUE)]
    allele2_tab_logr[,y_end := max(logr_dat$allele1_logr_no_dup, na.rm = TRUE)]
    
    allele2_logr_plot <- ggplot() +
      geom_rect(data  = allele2_tab_logr,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_line(data = logr_dat[allele_name == allele2_name],
                aes(allele1_snp_pos,
                    allele1_logr_no_dup), colour  = col_dt[allele_name == allele2_name]$allele_col) +
      geom_point(data = logr_dat[allele_name == allele2_name],
                 aes(allele1_snp_pos,
                     allele1_logr_no_dup),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      xlab("Position in  allele") +
      ylab("logR") +
      ggtitle(allele2_name)
    
    # plot the coverage
    allele1_tumour_coverage <- coverage_table_at_snps[,c("allele1_snp_pos",
                                                         "allele1_tumour_depth_no_dup")]
    allele1_tumour_coverage[,allele_name := allele1_name]
    allele1_tumour_coverage[,sample_type := "Tumour"]
    
    allele1_gl_coverage <- coverage_table_at_snps[,c("allele1_snp_pos",
                                                     "allele1_gl_depth_no_dup")]
    allele1_gl_coverage[,allele_name := allele1_name]
    allele1_gl_coverage[,sample_type := "Germline"]
    allele1_gl_coverage[,allele1_gl_depth_no_dup := allele1_gl_depth_no_dup * 1/mult_factor]
    
    allele2_tumour_coverage <- coverage_table_at_snps[,c("allele1_snp_pos",
                                                         "allele2_tumour_depth_no_dup")]
    allele2_tumour_coverage[,allele_name := allele2_name]
    allele2_tumour_coverage[,sample_type := "Tumour"]
    
    allele2_gl_coverage <- coverage_table_at_snps[,c("allele2_snp_pos",
                                                     "allele2_gl_depth_no_dup")]
    allele2_gl_coverage[,allele_name := allele2_name]
    allele2_gl_coverage[,sample_type := "Germline"]
    allele2_gl_coverage[,allele2_gl_depth_no_dup := allele2_gl_depth_no_dup * 1/mult_factor]
    
    coverage_tab <- rbindlist(list(allele1_tumour_coverage,
                                   allele1_gl_coverage,
                                   allele2_tumour_coverage,
                                   allele2_gl_coverage), use.names = FALSE)
    
    coverage_tab[allele_name == allele1_name & sample_type == "Tumour", colour := "#377EB8"]
    coverage_tab[allele_name == allele1_name & sample_type == "Germline", colour := "#A6CEE3"]
    coverage_tab[allele_name == allele2_name & sample_type == "Tumour", colour := "#E41A1C"]
    coverage_tab[allele_name == allele2_name & sample_type == "Germline", colour := "#FB9A99"]
    
    allele1_tab_cov <- gtf[V1 == allele1_name & V3 == "exon"]
    allele1_tab_cov[,y_start := min(coverage_tab$allele1_tumour_depth_no_dup, na.rm = TRUE)]
    allele1_tab_cov[,y_end := max(coverage_tab$allele1_tumour_depth_no_dup, na.rm = TRUE)]
    
    allele1_coverage_plot <- ggplot() +
      geom_rect(data  = allele1_tab_cov,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_line(data = coverage_tab[allele_name == allele1_name],
                aes(allele1_snp_pos,
                    allele1_tumour_depth_no_dup, colour = sample_type)) +
      geom_point(data = coverage_tab[allele_name == allele1_name],
                 aes(allele1_snp_pos,
                     allele1_tumour_depth_no_dup),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      scale_colour_manual(breaks =  c("Tumour", "Germline"),
                          values = c(unique(coverage_tab[allele_name == allele1_name & sample_type == "Tumour"]$colour), 
                                     unique(coverage_tab[allele_name == allele1_name & sample_type == "Germline"]$colour))) +
      xlab("Position in  allele") +
      ylab("Coverage") +
      ggtitle(allele1_name)
    
    allele2_tab_cov <- gtf[V1 == allele2_name & V3 == "exon"]
    allele2_tab_cov[,y_start := min(coverage_tab$allele1_tumour_depth_no_dup, na.rm = TRUE)]
    allele2_tab_cov[,y_end := max(coverage_tab$allele1_tumour_depth_no_dup, na.rm = TRUE)]
    
    allele2_coverage_plot <- ggplot() +
      geom_rect(data  = allele2_tab_cov,
                aes(xmin = V4,
                    xmax = V5,
                    ymin = y_start,
                    ymax = y_end,
                    fill = V3),
                alpha = 0.2) +
      geom_line(data = coverage_tab[allele_name == allele2_name],
                aes(allele1_snp_pos,
                    allele1_tumour_depth_no_dup, colour = sample_type)) +
      geom_point(data = coverage_tab[allele_name == allele2_name],
                 aes(allele1_snp_pos,
                     allele1_tumour_depth_no_dup),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank()) +
      scale_fill_manual(breaks =  c("exon"),
                        values = c("#D95F02")) +
      scale_colour_manual(breaks =  c("Tumour", "Germline"),
                          values = c(unique(coverage_tab[allele_name == allele2_name & sample_type == "Tumour"]$colour), 
                                     unique(coverage_tab[allele_name == allele2_name & sample_type == "Germline"]$colour))) +
      xlab("Position in  allele") +
      ylab("Coverage") +
      ggtitle(allele2_name)
  }else{
    cn_boxplot_plot <- ggplot()
    allele1_cn_plot <- ggplot()
    allele2_cn_plot <- ggplot()
    allele1_logr_plot <- ggplot()
    allele2_logr_plot <- ggplot()
    allele1_coverage_plot <- ggplot()
    allele2_coverage_plot <- ggplot()
  }
  
  return(list(cn1 = median_cn1,
              cn1_lower = cn1_lower,
              cn1_upper = cn1_upper,
              cn2 = median_cn2,
              cn2_lower = cn2_lower,
              cn2_upper = cn2_upper,
              cn1_binned = median_cn1_binned,
              cn1_binned_lower = cn1_binned_lower,
              cn1_binned_upper = cn1_binned_upper,
              cn2_binned  = median_cn2_binned,
              cn2_binned_lower = cn2_binned_lower,
              cn2_binned_upper = cn2_binned_upper,
              cn_boxplot_plot = cn_boxplot_plot,
              cn_bin_boxplot_plot =  cn_bin_boxplot_plot,
              allele1_cn_plot = allele1_cn_plot,
              allele2_cn_plot = allele2_cn_plot,
              allele1_logr_plot = allele1_logr_plot,
              allele2_logr_plot = allele2_logr_plot,
              allele1_coverage_plot = allele1_coverage_plot,
              allele2_coverage_plot = allele2_coverage_plot,
              n_bins = n_bins,
              n_snps = n_snps))
  
}

calculate_depth_aib  <- function(allele1_name,
                                 allele2_name,
                                 allele1_mismatch_bed_file,
                                 allele2_mismatch_bed_file,
                                 allele1_unique_coverage_file,
                                 allele2_unique_coverage_file,
                                 plot_title,
                                 make_plot){
  
  allele1_mismatch <- fread(allele1_mismatch_bed_file, sep = "\t")
  allele2_mismatch <- fread(allele2_mismatch_bed_file, sep = "\t")
  
  if(nrow(allele1_mismatch) == 0 | nrow(allele2_mismatch) == 0){
    return(list(paired_t_test = NA,
                paired_wilcoxon_test = NA,
                boxplot = ggplot(),
                allele1_median_dp = NA,
                allele2_median_dp = NA,
                allele1_n_snps_with_coverage = 0,
                allele2_n_snps_with_coverage = 0))
  }
  
  setnames(allele1_mismatch, c("CHROM", "chromStart", "chromEnd"))
  setnames(allele2_mismatch, c("CHROM", "chromStart", "chromEnd"))
  
  allele1_coverage <- fread(allele1_unique_coverage_file)
  allele2_coverage <- fread(allele2_unique_coverage_file)
  
  if(nrow(allele1_coverage) == 0 | nrow(allele2_coverage) == 0){
    # return(list(paired_t_test = NA,
    #             paired_wilcoxon_test = NA,
    #             boxplot = ggplot(),
    #             allele1_median_dp = NA,
    #             allele2_median_dp = NA,
    #             allele1_n_snps_with_coverage = nrow(allele1_coverage),
    #             allele2_n_snps_with_coverage = nrow(allele2_coverage)))
    stop("Why are there no lines?")
  }
  
  setnames(allele1_coverage, c("position", "allele1_depth"))
  setnames(allele2_coverage, c("position", "allele2_depth"))
  
  coverage_table_at_mismatch <- data.table(allele1_mismatch_pos = allele1_mismatch$chromEnd,
                                           allele2_mismatch_pos = allele2_mismatch$chromEnd)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele1_coverage,
                                      by.x = "allele1_mismatch_pos", by.y = "position", 
                                      all.x = TRUE)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele2_coverage,
                                      by.x = "allele2_mismatch_pos", by.y = "position", 
                                      all.x = TRUE)
  
  if(nrow(coverage_table_at_mismatch[is.na(allele1_depth) | is.na(allele2_depth) ]) > 0){
    stop("why are there missing lines?")
  }
  
  coverage_table_at_mismatch[,allele1_depth := as.numeric(allele1_depth)]
  coverage_table_at_mismatch[,allele2_depth := as.numeric(allele2_depth)]
  
  allele1_median_cov <- median(coverage_table_at_mismatch$allele1_depth, na.rm = TRUE)
  allele2_median_cov <- median(coverage_table_at_mismatch$allele2_depth, na.rm = TRUE)
  
  # if one of the alleles has duplicated mismatches, the depth for the other allele
  # will change, so we take the median
  coverage_table_at_mismatch[, allele1_depth_no_dup := median(allele1_depth), by = allele2_mismatch_pos]
  coverage_table_at_mismatch[, allele2_depth_no_dup := median(allele2_depth), by = allele1_mismatch_pos]
  
  # mark and remove duplicated values
  coverage_table_at_mismatch[, allele1_dup_id := 1:.N, allele1_mismatch_pos]
  coverage_table_at_mismatch[, allele2_dup_id := 1:.N, allele2_mismatch_pos]
  coverage_table_at_mismatch <- coverage_table_at_mismatch[allele1_dup_id == 1 & allele2_dup_id == 1]
  
  allele1_n_snps_with_coverage <- nrow(coverage_table_at_mismatch[allele1_depth_no_dup > 0])
  allele2_n_snps_with_coverage <- nrow(coverage_table_at_mismatch[allele2_depth_no_dup > 0])
  
  # get p values
  paired_t_test <- myTryCatch(t.test(coverage_table_at_mismatch$allele1_depth_no_dup,
                                     coverage_table_at_mismatch$allele2_depth_no_dup,
                                     paired=TRUE))
  
  paired_wilcoxon_test <- myTryCatch(wilcox.test(x = coverage_table_at_mismatch$allele1_depth_no_dup,
                                                 y = coverage_table_at_mismatch$allele2_depth_no_dup,
                                                 paired=TRUE))
  if(is.null(paired_t_test$value)){
    rna_aib_t_test <- NA
  }else{
    rna_aib_t_test <- paired_t_test$value$p.value
  }
  
  if(is.null(paired_wilcoxon_test$value)){
    rna_aib_wilcox_test <- NA
  }else{
    rna_aib_wilcox_test <- paired_wilcoxon_test$value$p.value
  }
  
  # plot the depth
  if(make_plot){
    allele1_depth_plot <- coverage_table_at_mismatch[,c("allele1_mismatch_pos",
                                                        "allele1_depth")]
    allele1_depth_plot[,allele_name := allele1_name]
    allele1_depth_plot[,allele_col := "#377EB8"]
    
    allele2_depth_plot <- coverage_table_at_mismatch[,c("allele1_mismatch_pos",
                                                        "allele2_depth")]
    allele2_depth_plot[,allele_name := allele2_name]
    allele2_depth_plot[,allele_col := "#E41A1C"]
    
    depth_plot_dat <- rbindlist(list(allele1_depth_plot, 
                                     allele2_depth_plot),
                                use.names = FALSE)
    
    col_dt <- unique(depth_plot_dat[,c("allele_name", "allele_col")])
    
    rna_aib_boxplot_plot <- ggplot(depth_plot_dat, 
                                   aes(allele_name, 
                                       allele1_depth,
                                       colour = allele_name)) +
      geom_boxplot() +
      geom_quasirandom() +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.position = "none") +
      ylab("Depth") +
      ggtitle(plot_title, 
              subtitle = paste0("t test p value=", 
                                formatC(paired_t_test$value$p.value, format = "e", digits = 2),
                                "\nwilcoxon p value=", 
                                formatC(paired_wilcoxon_test$value$p.value, format = "e", digits = 2))) +
      scale_colour_manual(breaks = c(col_dt$allele_name),
                          values = c(col_dt$allele_col))
  }else{
    rna_aib_boxplot_plot <- ggplot()
  }
  
  return(list(paired_t_test = rna_aib_t_test,
              paired_wilcoxon_test = rna_aib_wilcox_test,
              boxplot = rna_aib_boxplot_plot,
              allele1_median_dp = allele1_median_cov,
              allele2_median_dp = allele2_median_cov,
              allele1_n_snps_with_coverage = allele1_n_snps_with_coverage,
              allele2_n_snps_with_coverage = allele2_n_snps_with_coverage))
  
}

calculate_exp_dp <- function(allele1_name,
                             allele2_name,
                             allele1_gl_coverage_file,
                             allele2_gl_coverage_file,
                             allele1_snp_bed_file,
                             allele2_snp_bed_file,
                             purity,
                             gl_library_size,
                             tumour_library_size){
  
  allele1_snp <- fread(allele1_snp_bed_file, sep = "\t")
  allele2_snp <- fread(allele2_snp_bed_file, sep = "\t")
  
  if(nrow(allele1_snp) == 0 | nrow(allele2_snp) == 0){
    return(list(allele1_exp_dp = as.numeric(NA),
                allele2_exp_dp = as.numeric(NA),
                exp_dp_boxplot_plot = ggplot()))
  }
  
  setnames(allele1_snp, c("chr", "start", "end"))
  setnames(allele2_snp, c("chr", "start", "end"))
  
  allele1_gl_coverage <- fread(allele1_gl_coverage_file)
  setnames(allele1_gl_coverage, c("position", "allele1_gl_depth"))
  
  allele2_gl_coverage <- fread(allele2_gl_coverage_file)
  setnames(allele2_gl_coverage, c("position", "allele2_gl_depth"))
  
  mult_factor <- tumour_library_size/gl_library_size
  
  # create a table with the coverage at the snp positions
  coverage_table_at_snps <- data.table(allele1_snp_pos = allele1_snp[,end],
                                       allele2_snp_pos = allele2_snp[,end])
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele1_gl_coverage,
                                  by.x = "allele1_snp_pos", by.y = "position", 
                                  all.x = TRUE)
  
  coverage_table_at_snps <- merge(coverage_table_at_snps, allele2_gl_coverage,
                                  by.x = "allele2_snp_pos", by.y = "position", all.x = TRUE)
  
  coverage_table_at_snps[,allele1_exp_depth := allele1_gl_depth * mult_factor * purity]
  coverage_table_at_snps[,allele2_exp_depth := allele2_gl_depth * mult_factor * purity]
  
  # plot the expected dp
  allele1_dp <- coverage_table_at_snps[,c("allele1_snp_pos",
                                          "allele1_exp_depth")]
  allele1_dp[,allele_name := allele1_name]
  allele1_dp[,allele_col := "#377EB8"]
  
  allele2_dp <- coverage_table_at_snps[,c("allele1_snp_pos",
                                          "allele2_exp_depth")]
  
  allele2_dp[,allele_name := allele2_name]
  allele2_dp[,allele_col := "#E41A1C"]
  
  dp_dat <- rbindlist(list(allele1_dp,
                           allele2_dp),
                      use.names = FALSE)
  
  col_dt <- unique(dp_dat[,c("allele_name", "allele_col")])
  
  dp_boxplot_plot <- ggplot(dp_dat,
                            aes(allele_name,
                                allele1_exp_depth,
                                colour = allele_name)) +
    geom_boxplot() +
    geom_quasirandom() +
    theme_bw() +
    theme(text = element_text(size= 20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
          legend.position = "none") +
    ylab("Expected depth") +
    scale_colour_manual(breaks = c(col_dt$allele_name),
                        values = c(col_dt$allele_col))
  
  return(list(allele1_exp_dp = median(coverage_table_at_snps[,allele1_exp_depth], na.rm = TRUE),
              allele2_exp_dp = median(coverage_table_at_snps[,allele2_exp_depth], na.rm = TRUE),
              exp_dp_boxplot_plot = dp_boxplot_plot))
  
}

calculate_tumour_normal_comparison <- function(allele_name,
                                               gl_library_size,
                                               tumour_library_size,
                                               allele_mismatch_bed_file,
                                               allele_gl_coverage_file,
                                               allele_tumour_coverage_file,
                                               make_plot){
  

  mult_factor <- gl_library_size/tumour_library_size
  
  allele_mismatch <- fread(allele_mismatch_bed_file, sep = "\t")
  
  if(nrow(allele_mismatch) == 0){
    return(list(paired_t_test = NA,
                paired_wilcoxon_test = NA,
                boxplot_plot = ggplot(),
                median_tumour_rna_dp = NA,
                median_normal_rna_dp = NA,
                normal_n_snps_with_coverage = 0,
                tumour_n_snps_with_coverage = 0,
                coverage_plot = ggplot()))
  }
  
  setnames(allele_mismatch, c("CHROM", "chromStart", "chromEnd"))
  
  allele_normal_coverage <- fread(allele_gl_coverage_file)
  setnames(allele_normal_coverage, c("position", "allele_normal_depth"))
  
  allele_tumour_coverage <- fread(allele_tumour_coverage_file)
  setnames(allele_tumour_coverage, c("position", "allele_tumour_depth"))
  
  coverage_table_at_mismatch <- data.table(allele_mismatch_pos = allele_mismatch$chromEnd)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele_normal_coverage,
                                      by.x = "allele_mismatch_pos", by.y = "position", 
                                      all.x = TRUE)
  
  coverage_table_at_mismatch <- merge(coverage_table_at_mismatch, allele_tumour_coverage,
                                      by.x = "allele_mismatch_pos", by.y = "position", 
                                      all.x = TRUE)
  
  # normalise for differences in depth
  coverage_table_at_mismatch[,allele_tumour_normalised_depth := mult_factor * allele_tumour_depth]
  
  normal_n_snps_with_coverage <- nrow(coverage_table_at_mismatch[allele_normal_depth > 0])
  tumour_n_snps_with_coverage <- nrow(coverage_table_at_mismatch[allele_tumour_normalised_depth > 0])
  
  # get p values
  paired_t_test <- myTryCatch(t.test(coverage_table_at_mismatch$allele_normal_depth,
                                     coverage_table_at_mismatch$allele_tumour_normalised_depth,
                                     paired=TRUE))
  
  paired_wilcoxon_test <- myTryCatch(wilcox.test(coverage_table_at_mismatch$allele_normal_depth,
                                                 coverage_table_at_mismatch$allele_tumour_normalised_depth,
                                                 paired=TRUE))
  
  if(is.null(paired_t_test$value)){
    rna_tumour_normal_t_test <- NA
  }else{
    rna_tumour_normal_t_test <- paired_t_test$value$p.value
  }
  
  if(is.null(paired_wilcoxon_test$value)){
    rna_tumour_normal_wilcox_test <- NA
  }else{
    rna_tumour_normal_wilcox_test <- paired_wilcoxon_test$value$p.value
  }
  
  # plot the depth
  tumour_depth_plot <- coverage_table_at_mismatch[,c("allele_mismatch_pos",
                                                     "allele_tumour_normalised_depth")]
  
  
  normal_depth_plot <- coverage_table_at_mismatch[,c("allele_mismatch_pos",
                                                     "allele_normal_depth")]
  tumour_depth_plot[,sample_type := "Tumour"]
  tumour_depth_plot[,sample_col := "#E31A1C"]
  normal_depth_plot[,sample_type := "Normal"]
  normal_depth_plot[,sample_col := "#FB9A99"]
  
  depth_plot_dat <- rbindlist(list(tumour_depth_plot, 
                                   normal_depth_plot),
                              use.names = FALSE)
  
  col_dt <- unique(depth_plot_dat[,c("sample_type", "sample_col")])
  
  if(make_plot){
    rna_depth_boxplot_plot <- ggplot(depth_plot_dat, 
                                     aes(sample_type, 
                                         allele_tumour_normalised_depth,
                                         colour = sample_type)) +
      geom_boxplot() +
      geom_quasirandom() +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.position = "none") +
      ylab("RNA depth") +
      ggtitle(paste0("RNA tumour normal comparison:\n",allele_name), 
              subtitle = paste0("t test p value=", 
                                formatC(paired_t_test$value$p.value, format = "e", digits = 2),
                                "\nwilcoxon p value=", 
                                formatC(paired_wilcoxon_test$value$p.value, format = "e", digits = 2))) +
      scale_colour_manual(breaks = c(col_dt$sample_type),
                          values = c(col_dt$sample_col))
    
    coverage_plot <- ggplot() +
      geom_line(data = depth_plot_dat, 
                aes(allele_mismatch_pos, 
                    allele_tumour_normalised_depth,
                    colour = sample_type)) +
      geom_point(data = depth_plot_dat, 
                 aes(allele_mismatch_pos, 
                     allele_tumour_normalised_depth),
                 size = 0.5) +
      theme_bw() +
      theme(text = element_text(size= 20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust  = 0.5),
            legend.title = element_blank(),
            legend.position = "bottom") +
      scale_fill_manual(breaks =  paste0("Exon ", 1:8),
                        values = c("#1B9E77", "#D95F02", "#7570B3",
                                   "#E7298A", "#66A61E", "#E6AB02", 
                                   "#A6761D", "#666666")) +
      xlab("Position in transcript") +
      ylab("Coverage") +
      scale_colour_manual(breaks = c(col_dt$sample_type),
                          values = c(col_dt$sample_col)) +
      ggtitle(allele_name)
  }else{
    rna_depth_boxplot_plot  <- ggplot()
    coverage_plot  <- ggplot()
  }
  
  return(list(paired_t_test = rna_tumour_normal_t_test,
              paired_wilcoxon_test = rna_tumour_normal_wilcox_test,
              boxplot_plot = rna_depth_boxplot_plot,
              median_tumour_rna_dp = median(coverage_table_at_mismatch$allele_tumour_normalised_depth, na.rm = TRUE),
              median_normal_rna_dp = median(coverage_table_at_mismatch$allele_normal_depth, na.rm = TRUE),
              normal_n_snps_with_coverage = normal_n_snps_with_coverage,
              tumour_n_snps_with_coverage = tumour_n_snps_with_coverage,
              coverage_plot = coverage_plot))
}