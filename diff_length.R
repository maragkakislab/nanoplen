#!/usr/bin/env Rscript

# Debug test:
# dev/diff_length.R -d scratch/plen_test_data.tab -b "control" -l F -m scratch/plen_test_metadata.tab

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-d","--data_path"),
                help="Path to data file, long format, columns: library id, gene/transcript id, length"),
    make_option(c("-m","--metadata_path"),
                help="Path to metadata file, columns: library id, condition, [additional columns]"),
    make_option(c("-t","--test"), default = "t",
                help="Statistical test to use (t:t-test, m:linear mixed model, ks:kolmogorov-smirov, w:wilcoxon) [default %default]"),
    make_option(c("-b","--baseline"), default = "Control",
                help="String to specify baseline category [default %default]"),
    make_option(c("-l","--logscale"), default = "TRUE",
                help="Convert length to log2 scale [default %default]"),
    make_option(c("-o","--ofile"), default = "stdout",
                help="Path to output file [default %default]")
)
#FIXME: option for custom model in the future
opt <- parse_args(OptionParser(option_list = option_list))

vars = c("data_path", "metadata_path","test","baseline","logscale","ofile")
for (i in 1:length(vars)) {
    assign(vars[i],opt[[vars[i]]])
}
delim <- "\t"  #Because this script is directly after nanoplen, we can control the output
 
if (!(test %in% c("t","w","ks","m"))) {
    stop(sprintf("Unsupported test: %s. Accepted options: t, m, w, ks", test))
}

if (test %in% c("ks","w")) {
    stop(sprintf("Test %s not available yet, coming soon!",test))
}

# Input length data file and metadata file
data_file <- read.delim(data_path, header = TRUE, sep = delim)
metadata <- read.delim(metadata_path, header = TRUE, sep = delim)

#FIXME: Can omit this step if previous output names columns consistently
colnames(data_file) = c("lib_id","name","length")
colnames(metadata)[1:2] = c("lib_id", "condition")

# Remove rows with reported length 0. They should not be there anyway.
data_file = data_file[data_file$length > 0,]

# Add condition column from metadata
data_file$condition = sapply(data_file$lib_id, function(x) {metadata$condition[metadata$lib_id == x]})

# Relevel data_file$condition to use baseline string
data_file = within(data_file, condition <- relevel(condition, ref = baseline))


diff_length_single = function(data_file_sub, test, logscale = TRUE) {
    out = c(NA,NA)
    if (logscale) {
        data_file_sub$length = log2(data_file_sub$length)
        est_head = "log2FC"
    } else {
        est_head = "meandiff"
    }
    tryCatch({  
        if (test == "t") {
            res = lm(length~condition, data = data_file_sub)
            out = summary(res)$coefficients[2,c(1,4)]
        } else if (test == "m") {
            res = lme4::lmer(length~condition + (1 | lib_id), data = data_file_sub)
            out = summary(res)$coefficients[2,c(1,3)]
            out[2] = 2*pt(abs(out[2]), df=nrow(data_file_sub)-2,lower.tail = FALSE)
        }} 
    ,
        error = {function(e) {warning(
        #FIXME: Was supposed to also show which gene/transcript but cannot extract with current algorithm
            sprintf("Error: NAs given"))
        }}
    )
    names(out) = c(est_head, "pvalue")
    
    return(out)
}

diff_length = function(data_file, test) {
    data_file_byname = split(data_file, data_file$name)
    
    # Loops over all subsets split by name
    out = lapply(data_file_byname, function(d) {diff_length_single(d, test)})
    out = data.frame(do.call(rbind, out))
    out$qvalue = p.adjust(out$pvalue,method = "BH")
    return(out)
}

outres = diff_length(data_file, test)
outres = cbind(rownames(outres),outres)
colnames(outres)[1] = "name"

if (ofile == "stdout") {
    write.table(outres,file=stdout(),sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
} else {
    write.table(outres, file = ofile,sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
}
