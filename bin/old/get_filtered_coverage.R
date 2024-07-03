suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--coverage_path',  nargs=1,
                    help='Path to coverage file',
                    required=TRUE)
parser$add_argument('--filtered_positions_path',  nargs=1,
                    help='Path to bed containing positions in allele that pass minimum coverage in the germline',
                    required=TRUE)
parser$add_argument('--out_csv',  nargs=1,
                    help='Name for output coverage file',
                    required=TRUE)

args <- parser$parse_args()

coverage_path <- args$coverage_path
filtered_positions_path <- args$filtered_positions_path
out_csv <- args$out_csv

coverage_dt <- fread(coverage_path)
filtered_positions <- fread(filtered_positions_path, sep = "\t")

if(nrow(filtered_positions) > 0){
  
  setnames(filtered_positions, c("chr", "start", "end"))
  filtered_coverage_dt <- coverage_dt[pos %in% filtered_positions[,end]]
  
}else{
  filtered_coverage_dt <- data.table(chr = character(),
                                     start = numeric(),
                                     end = numeric())
}

fwrite(filtered_coverage_dt,	
       file = out_csv)


