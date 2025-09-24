#!/usr/bin/env Rscript

# This script is a bash wrapper for the differential length analysis functions
# 

suppressPackageStartupMessages(library(optparse))
library(nanoplen)

option_list <- list(
    make_option(c("-d","--data_path"),
                help="Path to data file, long format, columns: library id, gene/transcript id, length"),
    make_option(c("-m","--metadata_path"),
                help="Path to metadata file, columns: library id, condition, [additional columns]"),
    make_option(c("-t","--test"), default = "t",
                help="Statistical test to use (t:t-test, m:linear mixed model, w:wilcoxon, s:spline) [default %default]"),
    make_option(c("-c","--condition"), default = NULL,
                help="Condition variable to test on [default uses second metadata column]"),
    make_option(c("-b","--baseline"), default = NULL,
                help="String to specify baseline category"),
    make_option(c("--min_filter"), default = 30,
                help="Minimum number of reads to be included in analysis [default %default]"),
    make_option(c("-l","--logscale"), action = "store_true",  default=FALSE,
                help="Convert length to log2 scale (TRUE/FALSE) [default %default]"),
    make_option(c("-p","--params"), default = NULL,
                help="Extra parameters to use for using linear regression methods, separated by +, no spaces (example: time+age+age*time) [default %default]"),
    make_option(c("-n","--norm"), action = "store_true",  default=FALSE,
                help="Normalize data by normalization group. Requires norm_group column in metadata where replicates have the same value with each other [default %default]"),
    make_option(c("-f","--filter_file"), default=NULL,
                help="A file with contigs that you want to pre-filter for [default %default]"),
    make_option(c("-o","--ofile"), default = "stdout",
                help="Path to output file [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

delim <- "\t"  #Because this script is directly after nanoplen, we can control the output

if (!(opt$test %in% c("t","w","m","s"))) {
    stop(sprintf("Unsupported test: %s. Accepted options: t, m, w, s", opt$test))
}

# Input length data file and metadata file
data_file <- read.delim(opt$data_path, header = TRUE, sep = delim)
metadata <- read.delim(opt$metadata_path, header = TRUE, sep = delim)

if (!is.null(opt$filter_file)) {
  filter_file = read.table(opt$filter_file, header = F, sep = "\t", stringsAsFactors = F)
  data_file = data_file[data_file[,2] %in% opt$filter_file[,1],]
}

# Global variable. Low priority
has_warning <<- FALSE
wfile = ifelse(opt$ofile == "stdout", "warnings.txt", sprintf("%s_warnings.txt", opt$ofile))
ww <- file(wfile, open = "wt")
sink(ww, type = "message")

condition = opt$condition
if (!is.null(condition)) {
    if (!(condition %in% colnames(metadata))) {
        stop("Condition variable not in metadata columns!")
    } else {
        #Puts condition variable as second 
        #This assumes the sample variable is always first in the metadata
        ind = which(condition==colnames(metadata))
        metadata = metadata[,c(1,ind,(2:ncol(metadata))[-(ind-1)])]
    }
}

outres = nanoplen(data_file = data_file,
                  metadata = metadata,
                  test = opt$test,
                  baseline = opt$baseline,
                  min_filter = opt$min_filter,
                  logscale = opt$logscale,
                  params = opt$params,
                  norm = opt$norm)



if (opt$ofile == "stdout") {
    write.table(outres, file=stdout(), sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
} else {
    write.table(outres, file = opt$ofile, sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
}

sink(type="message")
close(ww)
if (has_warning) {
    message(sprintf("Some tests errored, logged in %s", wfile))
}
