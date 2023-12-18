suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--coverage_path',  nargs=1,
                    help='Path to coverage file',
                    required=TRUE)
parser$add_argument('--out_bed_name',  nargs=1,
                    help='Path to output bed',
                    required=TRUE)
parser$add_argument('--min_depth',  nargs=1,
                    help='Minimum depth for a position to pass',
                    required=TRUE)
parser$add_argument('--allele',  nargs=1,
                    help='Allele',
                    required=TRUE)

args <- parser$parse_args()

coverage_path <- args$coverage_path
out_bed_name <- args$out_bed_name
min_depth <- as.numeric(args$min_depth)
allele <- args$allele

coverage_dt <- fread(coverage_path)

filtered_coverage_dt <- coverage_dt[depth > min_depth] 
filtered_coverage_dt[, chr := allele]

filtered_coverage_bed <- filtered_coverage_dt[,c("chr", "pos")]
setnames(filtered_coverage_bed, c("chr", "end"))
filtered_coverage_bed[, start := end - 1 ]
setcolorder(filtered_coverage_bed, c("chr", "start", "end"))

fwrite(filtered_coverage_bed,	
       file = out_bed_name, 	
       sep="\t", col.names = FALSE)
