
xml_to_dt <- function(allele_xml){
  
  allele_feature_table <- data.table()
  sequence <- as.character(allele_xml$allele$sequence$nucsequence)
  
  # check that the sequence and the total length are the same
  
  for(feature_idx in 1:length(allele_xml$allele$sequence)){
    
    if(names(allele_xml$allele$sequence[feature_idx]) != "feature"){
      next
    }
    
    if(allele_xml$allele$sequence[feature_idx]$feature$.attrs["name"]=="Translation"){
      next
    }
    
    start <- as.character(allele_xml$allele$sequence[feature_idx]$feature$SequenceCoordinates["start"])
    end <- as.character(allele_xml$allele$sequence[feature_idx]$feature$SequenceCoordinates["end"])
    id <- as.character(allele_xml$allele$sequence[feature_idx]$feature$.attrs["id"])
    order <- as.character(allele_xml$allele$sequence[feature_idx]$feature$.attrs["order"])
    feature_type <- as.character(allele_xml$allele$sequence[feature_idx]$feature$.attrs["featuretype"])
    name <- as.character(allele_xml$allele$sequence[feature_idx]$feature$.attrs["name"])
    
    feature_sequence = substr(sequence, as.numeric(start), as.numeric(end))        
    
    tmp_dt <- data.table(feature_type = feature_type,
                         feature_name = name,
                         feature_start = start,
                         feature_end = end,
                         feature_id = id,
                         feature_order = order,
                         feature_sequence = feature_sequence)
    
    allele_feature_table <- rbindlist(list(allele_feature_table, tmp_dt), use.names = TRUE)
    
  }
  return(allele_feature_table)
}

get_dist_matrix <- function(msf_dir, gene_to_run){
  
  gene_to_run <- gsub("HLA-", "", gene_to_run)
  
  if(gene_to_run %in% c("DRB3", "DRB4", "DRB5")){
    msf_path <- paste0(msf_dir,"/DRB345_nuc.msf")
  }else{
    msf_path <- paste0(msf_dir,"/", gene_to_run, "_nuc.msf")
  }
  
  if(!file.exists(msf_path)){
    cat("msf path\n")
    cat(msf_path)
    stop("msf path doesn't exist")
  }
  
  msf <- read.alignment(msf_path, format = "msf", forceToLower = FALSE)
  # Get distance matrix
  dist_mat = as.matrix(dist.alignment(msf))
  # because DRB3, DRB4, DRB5 are all in the same msf file, we need to filter to
  # so the dist_mat only contains alleles from gene - otherwise we could be selecting
  # alleles from a different gene to be the similar allele
  dist_mat_rownames <- rownames(dist_mat)
  dist_mat_colnames <- colnames(dist_mat)
  dist_mat <- dist_mat[grepl(gene_to_run, dist_mat_rownames), grepl(gene_to_run, dist_mat_colnames)]
  
  # check all alleles are included in dist_mat
  
  return(dist_mat)
  
}

get_gen_seq_strand_specific <- function(gen_seq, strand){
  if(strand == "-"){
    allele_seq <- rev(gen_seq)
    strand_aware_dt <- data.table(base = allele_seq)
    strand_aware_dt[base == "G", new_base := "C"]
    strand_aware_dt[base == "C", new_base := "G"]
    strand_aware_dt[base == "A", new_base := "T"]
    strand_aware_dt[base == "T", new_base := "A"]
    
    strand_aware_gen_seq <- strand_aware_dt$new_base
    
  }else{
    strand_aware_gen_seq <- gen_seq
  }
  return(strand_aware_gen_seq)
}

get_gtf <- function(allele_dt, allele, gene, strand){
  
  gtf <- allele_dt[type %in% c("exon", "UTR", "intron")]
  gtf[,seqname := allele]
  gtf[,source := "imgt"]
  gtf[,score := "."]
  gtf[,strand := strand]
  gtf[,frame := "."]
  gtf[,similar_allele_save_name := paste0("hla_", (tolower(gsub('\\*|:', '_', similar_allele))))]
  gtf <- gtf[,c("seqname", "source", "type", "start", 
                "end", "score", "strand", "frame", "feature_seq_in_imgt", 
                "similar_allele_save_name", "min_dist")]
  gtf[, start := as.numeric(start)]
  gtf[, end := as.numeric(end)]
  
  # check that the number of introns and exons - HLA-B
  # if(nrow(gtf[type == "exon"]) != 7 & gene == "HLA-B"){
  #   stop("There should be 7 exons for HLA-B")
  # }
  # if(nrow(gtf[type == "intron"]) != 6 & gene == "HLA-B"){
  #   stop("There should be 6 introns for HLA-B")
  # }
  # 
  # # check that the number of introns and exons - HLA-A and HLA-C
  # if(nrow(gtf[type == "exon"]) != 8 & gene %in% c("HLA-C", "HLA-A")){
  #   stop("There should be 7 exons for HLA-B")
  # }
  # if(nrow(gtf[type == "intron"]) != 7 & gene %in% c("HLA-C", "HLA-A")){
  #   stop("There should be 6 introns for HLA-B")
  # }
  # if(!gene %in% c("HLA-B", "HLA-A", "HLA-C")){
  #   stop("Need to update")
  # }
  
  if(strand == "-"){
    
    coords <- data.table(old_coord = 1:max(gtf$end))
    coords[,new_coord  := rev(coords$old_coord)]
    
    strand_aware_gtf <- merge(gtf, coords, by.x = "start", by.y = "old_coord", all.x =TRUE)
    setnames(strand_aware_gtf, "new_coord","strand_aware_end")
    
    strand_aware_gtf <- merge(strand_aware_gtf, coords, by.x = "end", by.y = "old_coord", all.x =TRUE)
    setnames(strand_aware_gtf, "new_coord","strand_aware_start")
    
    strand_aware_gtf[,c("start", "end") := NULL]
    setnames(strand_aware_gtf, 
             c("strand_aware_start", "strand_aware_end"),
             c("feature_start", "feature_end"))
    
    # on reverse strand 3UTR will be first and 5UTR will be last
    strand_aware_gtf[,min_coord := min(feature_start), by = "type"]
    strand_aware_gtf[min_coord == feature_start & type == "UTR", type := "3UTR"]
    strand_aware_gtf[min_coord != feature_start & type == "UTR", type := "5UTR"]
    strand_aware_gtf[,min_coord := NULL]
    
    # feature name is counted from 5UTR
    # for reverse strand order by -pos
    strand_aware_gtf <- strand_aware_gtf[order(-feature_start)]
    strand_aware_gtf[, feature_number := 1:.N,  by = "type"]
    strand_aware_gtf[type == "exon", feature_name := paste0("exon_", feature_number)]
    strand_aware_gtf[type == "intron", feature_name := paste0("intron_", feature_number)]
    strand_aware_gtf[is.na(feature_name),feature_name := type]
    strand_aware_gtf[,feature_number := NULL]
    
  }else{
    strand_aware_gtf <- copy(gtf)
    setnames(strand_aware_gtf, 
             c("start", "end"),
             c("feature_start", "feature_end"))
    # on forward strand 5UTR will be first and 3UTR will be last
    strand_aware_gtf[,min_coord := min(feature_start), by = "type"]
    strand_aware_gtf[min_coord == feature_start & type == "UTR", type := "5UTR"]
    strand_aware_gtf[min_coord != feature_start & type == "UTR", type := "3UTR"]
    strand_aware_gtf[,min_coord := NULL]
    
    # exon name is counted from 5UTR
    # for forward strand order by pos
    strand_aware_gtf <- strand_aware_gtf[order(feature_start)]
    strand_aware_gtf[, feature_number := 1:.N,  by = "type"]
    strand_aware_gtf[type == "exon", feature_name := paste0("exon_", feature_number)]
    strand_aware_gtf[type == "intron", feature_name := paste0("intron_", feature_number)]
    strand_aware_gtf[is.na(feature_name),feature_name := type]
    strand_aware_gtf[,feature_number := NULL]
    
  }
  
  # add in the exon number which is just counted from pos
  strand_aware_gtf <- strand_aware_gtf[order(feature_start)]
  strand_aware_gtf[,exon_number := 1:.N, by = "type"]
  
  strand_aware_gtf[type == "exon" , 
                   other_values := paste0('gene_id "', gene,
                                          '"; transcript_id "', allele, 
                                          '"; exon_number "', exon_number,
                                          '"; feature_name "', feature_name,
                                          '"; exon_id "', allele, '.', exon_number,
                                          '"; feature_seq_in_imgt "', feature_seq_in_imgt,
                                          '"; similar_allele_save_name "', similar_allele_save_name,
                                          '"; min_dist "', min_dist,
                                          '"; gene_name "', gene, '";')]
  strand_aware_gtf[type != "exon" , 
                   other_values := paste0('gene_id "', gene,
                                          '"; transcript_id "', allele, 
                                          '"; feature_name "', feature_name,
                                          '"; exon_id "', allele, '.', exon_number,
                                          '"; feature_seq_in_imgt "', feature_seq_in_imgt,
                                          '"; similar_allele_save_name "', similar_allele_save_name,
                                          '"; min_dist "', min_dist,
                                          '"; gene_name "', gene, '";')]
  
  strand_aware_gtf <- strand_aware_gtf[,c("seqname", "source", "type",
                                          "feature_start", "feature_end",
                                          "score", "strand", "frame", "other_values")]
  
  return(strand_aware_gtf)
}

get_bed <- function(allele_dt, allele, gene, strand){
  
  allele_dt[,chr := allele]
  bed <- allele_dt[,c("chr", "start", "end", "type")]
  bed[, start := as.numeric(start)]
  bed[, end := as.numeric(end)]
  
  # check that the number of introns and exons - HLA-B
  if(nrow(bed[type == "exon"]) != 7 & gene == "HLA-B"){
    stop("There should be 7 exons for HLA-B")
  }
  if(nrow(bed[type == "intron"]) != 6 & gene == "HLA-B"){
    stop("There should be 6 introns for HLA-B")
  }
  
  # check that the number of introns and exons - HLA-A and HLA-C
  if(nrow(bed[type == "exon"]) != 8 & gene %in% c("HLA-C", "HLA-A")){
    stop("There should be 7 exons for HLA-B")
  }
  if(nrow(bed[type == "intron"]) != 7 & gene %in% c("HLA-C", "HLA-A")){
    stop("There should be 6 introns for HLA-B")
  }
  if(!gene %in% c("HLA-B", "HLA-A", "HLA-C")){
    stop("Need to update")
  }
  
  if(strand == "-"){
    
    coords <- data.table(old_coord = 1:max(bed$end))
    coords[,new_coord  := rev(coords$old_coord)]
    
    strand_aware_bed <- merge(bed, coords, by.x = "start", by.y = "old_coord", all.x =TRUE)
    setnames(strand_aware_bed, "new_coord","strand_aware_end")
    
    strand_aware_bed <- merge(strand_aware_bed, coords, by.x = "end", by.y = "old_coord", all.x =TRUE)
    setnames(strand_aware_bed, "new_coord","strand_aware_start")
    
    strand_aware_bed[,c("start", "end") := NULL]
    setnames(strand_aware_bed, 
             c("strand_aware_start", "strand_aware_end"),
             c("start", "end"))
    
    # on reverse strand 3UTR will be first and 5UTR will be last
    strand_aware_bed[,min_coord := min(start), by = "type"]
    strand_aware_bed[min_coord == start & type == "UTR", type := "3UTR"]
    strand_aware_bed[min_coord != start & type == "UTR", type := "5UTR"]
    strand_aware_bed[,min_coord := NULL]
    
    # feature name is counted from 5UTR
    # for reverse strand order by -pos
    strand_aware_bed <- strand_aware_bed[order(-start)]
    strand_aware_bed[, feature_number := 1:.N,  by = "type"]
    strand_aware_bed[type == "exon", feature_name := paste0("exon_", feature_number)]
    strand_aware_bed[type == "intron", feature_name := paste0("intron_", feature_number)]
    strand_aware_bed[is.na(feature_name),feature_name := type]
    strand_aware_bed[,feature_number := NULL]
    
  }else{
    strand_aware_gtf <- copy(gtf)
    setnames(strand_aware_gtf, 
             c("start", "end"),
             c("feature_start", "feature_end"))
    # on forward strand 5UTR will be first and 3UTR will be last
    strand_aware_gtf[,min_coord := min(feature_start), by = "type"]
    strand_aware_gtf[min_coord == feature_start & type == "UTR", type := "5UTR"]
    strand_aware_gtf[min_coord != feature_start & type == "UTR", type := "3UTR"]
    strand_aware_gtf[,min_coord := NULL]
    
    # exon name is counted from 5UTR
    # for forward strand order by pos
    strand_aware_gtf <- strand_aware_gtf[order(feature_start)]
    strand_aware_gtf[, feature_number := 1:.N,  by = "type"]
    strand_aware_gtf[type == "exon", feature_name := paste0("exon_", feature_number)]
    strand_aware_gtf[type == "intron", feature_name := paste0("intron_", feature_number)]
    strand_aware_gtf[is.na(feature_name),feature_name := type]
    strand_aware_gtf[,feature_number := NULL]
    
  }
  
  bed[, start := start - 1]
  
  # add in the exon number which is just counted from pos
  strand_aware_gtf <- strand_aware_gtf[order(feature_start)]
  strand_aware_gtf[,exon_number := 1:.N, by = "type"]
  
  strand_aware_gtf[type == "exon" , 
                   other_values := paste0('gene_id "', gene,
                                          '"; transcript_id "', allele, 
                                          '"; exon_number "', exon_number,
                                          '"; feature_name "', feature_name,
                                          '"; exon_id "', allele, '.', exon_number,
                                          '"; feature_seq_in_imgt "', feature_seq_in_imgt,
                                          '"; similar_allele_save_name "', similar_allele_save_name,
                                          '"; min_dist "', min_dist,
                                          '"; gene_name "', gene, '";')]
  strand_aware_gtf[type != "exon" , 
                   other_values := paste0('gene_id "', gene,
                                          '"; transcript_id "', allele, 
                                          '"; feature_name "', feature_name,
                                          '"; exon_id "', allele, '.', exon_number,
                                          '"; feature_seq_in_imgt "', feature_seq_in_imgt,
                                          '"; similar_allele_save_name "', similar_allele_save_name,
                                          '"; min_dist "', min_dist,
                                          '"; gene_name "', gene, '";')]
  
  strand_aware_gtf <- strand_aware_gtf[,c("seqname", "source", "type",
                                          "feature_start", "feature_end",
                                          "score", "strand", "frame", "other_values")]
  
  return(strand_aware_gtf)
}

get_allele_list <- function(allele_list_path){
  
  all_allele_list <- fread(allele_list_path)
  
  # check that alleles arent repeated:
  if(nrow(all_allele_list[,.N,by = "Allele"][N>1])){
    stop("Why are alleles repeated?")
  }
  
  all_allele_list[,gene := gsub("\\*.*", "", Allele)]
  all_allele_list[!gene %in% c("TAP1", "TAP2", "HFE", "MICA", "MICB"), gene := paste0("HLA-", gene)]
  all_allele_list[Partial == "Full" & Type == "gDNA", full_gen_seq_exists := TRUE]
  all_allele_list[Partial == "Full" & Type == "gDNA", full_nuc_seq_exists := TRUE]
  
  all_allele_list[Partial == "Full" & Type == "cDNA", full_gen_seq_exists := FALSE]
  all_allele_list[Partial == "Full" & Type == "cDNA", full_nuc_seq_exists := TRUE]
  
  all_allele_list[Partial == "Partial" & Type == "gDNA", full_gen_seq_exists := FALSE]
  all_allele_list[Partial == "Partial" & Type == "gDNA", full_nuc_seq_exists := FALSE]
  
  all_allele_list[Partial == "Partial" & Type == "cDNA", full_gen_seq_exists := FALSE]
  all_allele_list[Partial == "Partial" & Type == "cDNA", full_nuc_seq_exists := FALSE]
  
  return(all_allele_list)
}

get_missing_seq <- function(partial_allele, similar_allele, hla_dat_allele_features){
  
  if(length(similar_allele) != 1){
    stop("Similar allele should be length 1")
  }
  
  allele_dt <- hla_dat_allele_features[allele_name == partial_allele & type != "CDS"]
  similar_allele_dt <- hla_dat_allele_features[allele_name == similar_allele & type != "CDS"]
  
  if(nrow(allele_dt) == 0){
    stop("Allele not in hla_dat")
  }
  if(nrow(similar_allele_dt) == 0){
    stop("Similar allele not in hla_dat")
  }
  
  # check that the introns and exons have a number. 
  # this number comes from the hla.dat file
  # This is really important as we want to match "exon 1" with "exon 1" 
  # and the first exon that we have the sequence for might not be exon 1
  # note we dont have to worry about which direction the exon is counted from 
  # (which will differ depending on the strand) as both the partial allele 
  # and similar allele are from the same gene and so on the same strand
  if(nrow(allele_dt[type != "UTR" & is.na(number)]) > 0){
    stop("Missing feature numbers")
  }
  
  if(nrow(similar_allele_dt[type != "UTR" & is.na(number)]) > 0){
    stop("Missing feature numbers")
  }
  
  # add in feature numbers to UTRs
  similar_allele_dt <- similar_allele_dt[order(start)]
  similar_allele_dt[, utr_number := 1:.N, by = "type"]
  similar_allele_dt[type == "UTR", number := utr_number]
  similar_allele_dt[, utr_number := NULL]
  
  allele_dt <- allele_dt[order(start)]
  allele_dt[, utr_number := 1:.N, by = "type"]
  allele_dt[type == "UTR", number := utr_number]
  allele_dt[, utr_number := NULL]
  
  # add in the feature name to merge
  allele_dt[,feature_name := paste0(type, "_", number)]
  similar_allele_dt[,feature_name := paste0(type, "_", number)]
  
  # merge the two tables
  setnames(similar_allele_dt, paste0("similar_allele_", names(similar_allele_dt)))
  setnames(allele_dt, paste0("allele_", names(allele_dt)))
  
  both_alleles_dt <- merge(similar_allele_dt, allele_dt, 
                           by.x = "similar_allele_feature_name",
                           by.y = "allele_feature_name",
                           all = TRUE)
  
  # check that the similar allele isn't missing  anything
  if(nrow(both_alleles_dt[is.na(similar_allele_feature_seq)]) > 0){
    stop("Similar allele missing something")
  }
  
  # add the sequences that we do have to the similar allele
  both_alleles_dt[!is.na(allele_feature_seq), similar_allele_feature_seq := allele_feature_seq]
  
  # I need to update the start and end as the replaced
  # sequences might have different lengths - this will be used in the gtf later
  both_alleles_dt[,seq_length := nchar(similar_allele_feature_seq)]
  
  # I want to double check that the features in the similar allele don't overlap 
  # and don't contain any gaps
  
  bases_covered <- c()
  for(i in 1:nrow(both_alleles_dt)){
    bases_covered <- c(bases_covered, both_alleles_dt[i]$similar_allele_start:both_alleles_dt[i]$similar_allele_end)
  }
  
  # check bases aren't repeated
  if(length(unique(bases_covered)) != length(bases_covered)){
    stop("Features overlap")
  }
  if(any(!min(bases_covered):max(bases_covered) %in% bases_covered)){
    stop("Gaps between features")
  }
  
  # add in the start and end - we have to recalculate it as it might have changed when
  # we replaced feature sequences
  both_alleles_dt[, similar_allele_start := as.numeric(similar_allele_start)]
  both_alleles_dt <- both_alleles_dt[order(similar_allele_start)]
  both_alleles_dt[,feature_end := cumsum(both_alleles_dt$seq_length)]
  both_alleles_dt[,feature_start := feature_end - seq_length + 1]
  both_alleles_dt[!is.na(allele_feature_seq), feature_seq_in_imgt := TRUE]
  both_alleles_dt[is.na(allele_feature_seq), feature_seq_in_imgt := FALSE]
  
  allele_dt <- both_alleles_dt[,c("similar_allele_feature_name", "similar_allele_type",
                                  "similar_allele_number",
                                  "feature_start", "feature_end", "feature_seq_in_imgt",
                                  "similar_allele_feature_seq")]
  
  setnames(allele_dt,
           c("feature_name", "type",
             "number",
             "start", "end", "feature_seq_in_imgt",
             "feature_seq"))
  
  return(allele_dt)
  
}

get_cds_line <- function(ft_lines){
  
  # get codon start line
  cds_start_line_idx <- grep(".*[[:space:]]CDS[[:space:]].*", ft_lines)
  
  if(length(cds_start_line_idx) == 0){
    stop("No cds line?")
  }
  
  if(length(cds_start_line_idx) != 1){
    stop("There should only be one CDS line start")  
  }
  
  # loop over the cds lines until you get a "/"
  cds_lines <- c()
  for(cds_line_idx in cds_start_line_idx:length(ft_lines)){
    
    # check if it contains "/"
    stop_loop <- grepl('\\/', ft_lines[cds_line_idx])
    if(stop_loop){
      break
    }
    
    cds_lines <- c(cds_lines,
                   gsub(".*?([0-9].*[0-9]).*", "\\1", ft_lines[cds_line_idx]))
  }
  if(length(cds_lines) == 0){
    stop("Why no CDS lines?")
  }
  
  cds_line <- paste0(cds_lines, collapse = ",")
  return(cds_line)
  
}

get_translation_line <- function(ft_lines){
  
  # get codon start line
  translation_start_line_idx <- grep("translation=", ft_lines)
  
  if(length(translation_start_line_idx) > 1){
    stop("There should only be one translation line")
  }
  
  if(length(translation_start_line_idx) == 0){
    # stop("No translation line?")
    return(as.character(NA))
  }else{
    translation_start_line <- grep("translation=", ft_lines, value = TRUE)
    translation_start_line <- gsub('.*translation=\"', "", translation_start_line)
    translation_start_line_contains_end <- grepl('\"', translation_start_line)
    
    translation_next_lines <- c()
    if(!translation_start_line_contains_end){
      
      for(i in (translation_start_line_idx + 1):length(ft_lines)){
        
        translation_next_lines <- c(translation_next_lines, gsub("FT.*?([A-Z]+.*)", "\\1", ft_lines[i]))
        translation_next_line_contains_end <- grepl('\"', gsub("FT.*?([A-Z]+.*)", "\\1", ft_lines[i]))
        if(translation_next_line_contains_end){
          break
        }
      }
    }
    translation_line <- paste0(c(translation_start_line, translation_next_lines), collapse = "")
    translation_line <- gsub('\"', "", translation_line)
    return(translation_line)
  }
  
}

get_coords <-function(ft_lines){
  coord_dt <- data.table()
  for(i in 1:length(ft_lines)){
    
    if(grepl("[[:space:]]UTR[[:space:]]", ft_lines[i])){
      start_end <- gsub(".*UTR.*?([0-9+].*)", "\\1", ft_lines[i])
      coord_dt <- rbindlist(list(coord_dt, 
                                 data.table(start = gsub("([0-9]+)..[0-9]+", "\\1", start_end),
                                            end = gsub("[0-9]+..([0-9]+)", "\\1", start_end),
                                            type = "UTR",
                                            number = NA)))
    }
    
    if(grepl("[[:space:]]exon[[:space:]]", ft_lines[i])){
      start_end <- gsub(".*exon.*?([0-9+].*)", "\\1", ft_lines[i])
      number_line <- ft_lines[i + 1]
      number <- gsub('.*number=\"([0-9]+).*', "\\1", number_line)
      coord_dt <- rbindlist(list(coord_dt, 
                                 data.table(start = gsub("([0-9]+)..[0-9]+", "\\1", start_end),
                                            end = gsub("[0-9]+..([0-9]+)", "\\1", start_end),
                                            type = "exon",
                                            number = number)))
    }
    
    if(grepl("[[:space:]]intron[[:space:]]", ft_lines[i])){
      start_end <- gsub(".*intron.*?([0-9+].*)", "\\1", ft_lines[i])
      number_line <- ft_lines[i + 1]
      number <- gsub('.*number=\"([0-9]+).*', "\\1", number_line)
      coord_dt <- rbindlist(list(coord_dt, 
                                 data.table(start = gsub("([0-9]+)..[0-9]+", "\\1", start_end),
                                            end = gsub("[0-9]+..([0-9]+)", "\\1", start_end),
                                            type = "intron",
                                            number = number)))
    }
  }
  return(coord_dt)
}



