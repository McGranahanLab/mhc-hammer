suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--bam_path',  nargs=1,
                    help='Path to allele BAM',
                    required=TRUE)
parser$add_argument('--max_mismatch',  nargs=1,
                    help='Maximum number of mismatches',
                    required=TRUE)
parser$add_argument('--allele',  nargs=1,
                    help='Allele name',
                    required=TRUE)
parser$add_argument('--sample_name',  nargs=1,
                    help='Sample name',
                    required=TRUE)
parser$add_argument('--paired_end', nargs='?',
                    const=TRUE,
                    default=TRUE,
                    help=paste0("Specify whether the sample is paired-end. ",
                                "Include the argument (with or without a value) ",
                                "for paired-end data. If omitted, single-end ",
                                "data will be assumed."),
                    required=FALSE)

count.events <- function(BAMfile, n, pairedEnd){
  x              <- scanBam(BAMfile, index = BAMfile, param=ScanBamParam(what = scanBamWhat(), tag = 'NM'))
  readIDs        <- x[[1]][['qname']]
  cigar          <- x[[1]][['cigar']]
  editDistance   <- unlist(x[[1]][['tag']])
  insertionCount <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'I', x = unlist(strsplit(boop, split = '')))))} )
  deletionCount  <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'D', x = unlist(strsplit(boop, split = '')))))} )
  indelTotals    <- sapply(cigar, FUN = function(boop) {
    tmp <- unlist(strsplit(gsub("([0-9]+)","~\\1~",boop), "~" ))
    Is  <- grep(pattern = 'I', x = tmp)
    Ds  <- grep(pattern = 'D', x = tmp)
    total <- sum(as.numeric(tmp[(Is-1)])) + sum(as.numeric(tmp[Ds-1]))
    return(total)
  })
  misMatchCount <- editDistance - indelTotals
  eventCount <- misMatchCount + insertionCount + deletionCount
  names(eventCount) <- 1:length(eventCount)
  passed     <- eventCount[which(eventCount <= n)]
  y <- readIDs[as.numeric(names(passed))]
  if (pairedEnd) {
    y <- names(table(y)[which(table(y) == 2)])
  }
  return(y)
}

# parse args

args <- parser$parse_args()

bam_path <- args$bam_path
max_mismatch <- args$max_mismatch
allele <- args$allele
sample_name <- args$sample_name
paired_end <- args$paired_end

# Convert paired_end to logical if it's a string
if (is.character(paired_end)) {
  paired_end <- tolower(paired_end) == "true"
}

cat("bam_path=", bam_path, "\n")
cat("allele=", allele, "\n")
cat("sample_name=", sample_name, "\n")

if(!file.exists(bam_path)){
  stop("Bam doesn't exist: ", bam_path, "\n")
}

passed_reads <- count.events(bam_path, n = max_mismatch, pairedEnd = paired_end)
write.table(passed_reads, file = paste0(sample_name, ".", allele,".passed_reads.txt"),
                  sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
