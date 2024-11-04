library(argparse)
library(data.table)

parser <- ArgumentParser()

parser$add_argument('--sample_sheet_path',  nargs=1,
                    help='Path to sample sheet.',
                    required=TRUE)
parser$add_argument('--validated_sample_sheet_path', nargs=1,
                    help="Path to save sample sheet",
                    required=TRUE)
parser$add_argument('--run_hlahd', nargs='?',
                    const=TRUE,
                    default=TRUE,
                    help="If true, predict HLA alleles from WXS germline sample",
                    required=FALSE)

args <- parser$parse_args()

sample_sheet_path <- args$sample_sheet_path
validated_sample_sheet_path <- args$validated_sample_sheet_path
run_hlahd <- args$run_hlahd

# Read in the sample sheet
if(!file.exists(sample_sheet_path)){
  stop("File does not exist: ", sample_sheet_path)
}

sample_sheet <- fread(sample_sheet_path)

# check the columns that are required are present
required_col_names <- c("patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_name")
missing_cols <- any(!required_col_names %in% names(sample_sheet))

if(missing_cols){
  missing_col_names <- required_col_names[!required_col_names %in% names(sample_sheet)]
  msg <- paste0("Input sample sheet is missing the following columns:\n",
                paste0(missing_col_names, collapse = "\n"), "\n")
  stop(msg)
}

# check no duplicated samples
duplicated_samples <- nrow(sample_sheet[,.N,by = c("sample_name", "sequencing_type")][N>1]) > 0
if(duplicated_samples){
  duplicated_sample_names <- unique(sample_sheet[,.N,by = c("sample_name", "sequencing_type")][N>1]$sample_name)
  msg <- paste0("The following samples are duplicated (i.e. the sample name appears in more than one row with the same sequencing type):\n",
                paste0(duplicated_sample_names, collapse = "\n"), "\n")
  stop(msg)
}

# remove white space from patient and sample names
sample_sheet[,patient := gsub("\\s", "_", patient)]
sample_sheet[,sample_name := gsub("\\s", "_", sample_name)]
sample_sheet[,normal_sample_name := gsub("\\s", "_", normal_sample_name)]

# check bam ends with bam
bad_bams <- nrow(sample_sheet[!grepl("bam$", bam_path)]) > 0
if(bad_bams){
  bad_bam_paths <- sample_sheet[!grepl("bam$", bam_path)]$bam_path
  msg <- paste0("The following BAM files do not end in the suffix '.bam':\n",
                paste0(bad_bam_paths, collapse = "\n"), "\n")
  stop(msg)
}

# check bam file exists
sample_sheet[,bam_exists := file.exists(bam_path)]
missing_bams <- nrow(sample_sheet[bam_exists == FALSE]) > 0
if(missing_bams){
  missing_bam_paths <- sample_sheet[bam_exists == FALSE]$bam_path
  msg <- paste0("The following BAM files do not exist:\n",
                paste0(missing_bam_paths, collapse = "\n"), "\n")
  stop(msg)
}

# check sequencing type
sample_sheet[sequencing_type %in% c("wes", "wxs"), sequencing_type := "wxs"]
sample_sheet[sequencing_type %in% c("rnaseq", "rna"), sequencing_type := "rnaseq"]
bad_seq_type <- nrow(sample_sheet[!sequencing_type %in% c("wxs", "rnaseq")]) > 0
if(bad_seq_type){
  
  bad_seq_type_names <- unique(sample_sheet[!sequencing_type %in% c("wxs", "rnaseq")]$sequencing_type)
  msg <- paste0("The following entries for the sequencing_type column are not supported:\n",
                paste0(bad_seq_type_names, collapse = "\n"), "\n",
                "Please change to either 'wes', 'wxs', 'rnaseq', 'rna'.")
  stop(msg)
}

# check purity is between 0 and 1
bad_purity <- nrow(sample_sheet[purity > 1 | purity < 0]) > 0
if(bad_purity){
  
  bad_purity_samples <- unique(sample_sheet[purity > 1 | purity < 0]$sample_name)
  msg <- paste0("The following samples have a purity that does not lie between 0 and 1:\n",
                paste0(bad_purity_samples, collapse = "\n"), "\n")
  stop(msg)
}

# check every patient has a wxs germline sample
patients_with_wxs_gl <- unique(sample_sheet[sample_type == "normal" & sequencing_type == "wxs"]$patient)
missing_gl <- nrow(sample_sheet[!patient %in% patients_with_wxs_gl]) > 0
if(missing_gl){
  
  missing_gl_patients <- unique(sample_sheet[!patient %in% patients_with_wxs_gl]$patient)
  msg <- paste0("The following patients are missing a germline WXS sample in the sample sheet:\n",
                paste0(missing_gl_patients, collapse = "\n"), "\n")
  warning(msg)
  if (run_hlahd) {
    msg <- paste0("Cannot run HLAHD without a germline WXS, ",
                  "provide path to HLA alleles in input inventory file.")
    stop(msg)
  }
}

# check every matched normal has a bam path
matched_normals <- unique(sample_sheet[normal_sample_name != ""]$normal_sample_name)
missing_normals <- matched_normals[!matched_normals %in% sample_sheet$sample_name]
if(length(missing_normals) > 0){
  
  msg <- paste0("The following samples in the normal_sample_name column aren't in present in the sample_name column, and so do not have an associated bam_path:\n",
                paste0(missing_normals, collapse = "\n"), 
                "\nMake sure that every sample in normal_sample_name is also represented as a row in the sample sheet.\n")
  stop(msg)
}

fwrite(sample_sheet, validated_sample_sheet_path)
