suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--patient_id',  nargs=1,
                    help='Patient ID.',
                    required=TRUE)
parser$add_argument('--patient_hla_alleles', nargs=1,
                    help="Path to .csv file containing patient's HLA calls.",
                    required=TRUE)
parser$add_argument('--mhc_gtf',  nargs=1,
                    help='Path to HLA gtf',
                    required=TRUE)

args <- parser$parse_args()
patient_id <- args$patient_id
patient_hla_alleles_path <- args$patient_hla_alleles
mhc_gtf_path <- args$mhc_gtf

bed_path <- paste0(patient_id, ".bed")
exon_bed_path <- paste0(patient_id, ".exons.bed")

# get patient hla alleles
hlahd_tab <- fread(patient_hla_alleles_path, header = FALSE)
hlahd_tab <- hlahd_tab[V2 != "not typed" & V3 != "not typed"]
hlahd_tab <- hlahd_tab[V1 %in% c("A", "B", "C")]
patient_hla_alleles <- unique(c(hlahd_tab$V2, hlahd_tab$V3))

# read in gtf 
gtf <- fread(mhc_gtf_path)

# filter gtf for patient specific alleles
gtf <- gtf[seqname %in% patient_hla_alleles]
gtf[, feature_name := gsub('.*feature_name "(.*?)".*', "\\1", other_values)]

# create bed 
bed <- gtf[,c("seqname", "feature_start", "feature_end", "feature_name")]
bed[,feature_start := feature_start-1]

# write bed
fwrite(bed, file = bed_path, sep = "\t", col.names = FALSE, quote = FALSE)

# make a bed with just the exons
exon_bed <- gtf[type == "exon",
                c("seqname", "feature_start", "feature_end", "feature_name")]
exon_bed[, exon_number := 1:.N, by = "seqname"]
exon_bed[,width := feature_end - feature_start]

# update the start and end coordinates
for(seq in unique(exon_bed$seqname)){
  for(exon_num in sort(exon_bed[seqname == seq]$exon_number)){
    if(exon_num == 1){
      exon_bed[seqname == seq & exon_number == exon_num, new_start := feature_start]
      exon_bed[seqname == seq & exon_number == exon_num, new_end := new_start + width]
    }else{
      previous_exon_end <- exon_bed[seqname == seq & exon_number == (exon_num - 1)]$new_end
      exon_bed[seqname == seq & exon_number == exon_num, new_start := (previous_exon_end + 1)]
      exon_bed[seqname == seq & exon_number == exon_num, new_end := new_start + width]
    }
  }
}
exon_bed[, min_start := min(new_start), by = "seqname"]
exon_bed[, new_start := new_start - min_start + 1]
exon_bed[, new_end := new_end - min_start + 1]

exon_bed <- exon_bed[,c("seqname", "new_start", "new_end", "feature_name")]
setnames(exon_bed, c("seqname", "feature_start", "feature_end", "feature_name"))
exon_bed[,feature_start := feature_start-1]

# write bed
fwrite(exon_bed, file = exon_bed_path, sep = "\t", col.names = FALSE, quote = FALSE)


