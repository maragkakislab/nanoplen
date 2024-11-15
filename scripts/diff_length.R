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

vars = sapply(1:(length(option_list)), function(i) {substr(option_list[[i]]@long_flag,3,1000)})
for (i in 1:length(vars)) {
    assign(vars[i],opt[[vars[i]]])
}
delim <- "\t"  #Because this script is directly after nanoplen, we can control the output

if (!(test %in% c("t","w","m","s"))) {
    stop(sprintf("Unsupported test: %s. Accepted options: t, m, w, s", test))
}

# Input length data file and metadata file
data_file <- read.delim(data_path, header = TRUE, sep = delim)
metadata <- read.delim(metadata_path, header = TRUE, sep = delim)

if (!is.null(filter_file)) {
  filter_file = read.table(filter_file, header = T, sep = "\t", stringsAsFactors = F)
  data_file = data_file[data_file[,1] %in% filter_file[,1],]
}

# Global variable. Low priority
has_warning <<- FALSE
wfile = ifelse(ofile == "stdout", "warnings.txt", sprintf("%s_warnings.txt", ofile))
ww <- file(wfile, open = "wt")
sink(ww, type = "message")

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

outres = nanoplen(data_file,
                  metadata,
                  test,
                  baseline,
                  logscale,
                  params,
                  norm)



if (ofile == "stdout") {
    write.table(outres, file=stdout(), sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
} else {
    write.table(outres, file = ofile, sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
}

sink(type="message")
close(ww)
if (has_warning) {
    message(sprintf("Some tests errored, logged in %s", wfile))
}
