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
                help="Statistical test to use (t:t-test, m:linear mixed model, w:wilcoxon) [default %default]"),
    make_option(c("-b","--baseline"), default = "Control",
                help="String to specify baseline category [default %default]"),
    make_option(c("-l","--logscale"), action = "store_true",  default=FALSE,
                help="Convert length to log2 scale (TRUE/FALSE) [default %default]"),
    make_option(c("-p","--params"), default = NULL,
                help="Extra parameters to use for using linear regression methods, separated by +, no spaces (example: time+age+age*time) [default %default]"),
    make_option(c("-n","--norm"), action = "store_true",  default=FALSE,
                help="Normalize data by normalization group. Requires norm_group column in metadata where replicates have the same value with each other [default %default]"),
    make_option(c("-o","--ofile"), default = "stdout",
                help="Path to output file [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

vars = c("data_path", "metadata_path","test","baseline","logscale","params", "norm", "ofile")
for (i in 1:length(vars)) {
    assign(vars[i],opt[[vars[i]]])
}
delim <- "\t"  #Because this script is directly after nanoplen, we can control the output

if (!(test %in% c("t","w","m"))) {
    stop(sprintf("Unsupported test: %s. Accepted options: t, m, w", test))
}


# Input length data file and metadata file
data_file <- read.delim(data_path, header = TRUE, sep = delim)
metadata <- read.delim(metadata_path, header = TRUE, sep = delim)

# Global variable. Low priority
has_warning <<- FALSE
wfile = ifelse(ofile == "stdout", "warnings.txt", sprintf("%s_warnings.txt", ofile))
ww <- file(wfile, open = "wt")
sink(ww, type = "message")

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
