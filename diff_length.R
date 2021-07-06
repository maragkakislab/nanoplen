#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
    make_option(c("-p","--data"),
                help="Path to data file"),
    make_option(c("-m","--metadata"),
                help="Path to metadata file"),
    make_option(c("-t","--test"), default = "t",
                help="Statistical test to use [default %default]"),
    make_option(c("-b","--baseline"), default = "Control",
                help="String to specify baseline category [default %default]"),
    make_option(c("-l","--logscale"), default = "TRUE",
                help="Convert  [default %default]"),
    make_option(c("-o","--ofile"), default = "stdout()",
                help="Path to output file [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

vars = c("data_path", "metadata_path","test","baseline","logscale","ofile")
for (i in length(vars)) {
    assign(vars[i],opt[i])
}
delim <- "\t"  #Because this script is directly after nanoplen, we can control the output

if (!(test %in% c("t","w","ks","m"))) {
    stop(sprintf("Unsupported test: %s", test))
}

if (test %in% c("ks","w")) {
    stop(sprintf("Test %s not available yet, coming soon!",test))
}

#Input length data file and metadata file
data_file <- read.delim(data_path, header = TRUE, sep = delim)
metadata <- read.delim(metadata_path, header = TRUE, sep = delim)

#FIXME: Can omit this step if previous output names columns consistently
colnames(data_file) = c("lib_id","transcript","Y")

#Add condition column from metadata
data_file$condition = sapply(data_file[,1], function(x) {metadata[metadata[,1] == x,2]})

# Relevel data_file$condition to use baseline string
data_file = within(data_file, condition <- relevel(condition, ref = baseline))


diff_length_single = function(data_file_sub, test, logscale = TRUE) {
    #FIXME: tryCatch for potential errors
    out = c(NA,NA)
    if (logscale) {
        data_file_sub$Y = log2(data_file_sub$Y+1)
        est_head = "log2FC"
    } else {
        est_head = "meandiff"
    }
    tryCatch({  
        if (test == "t") {
            res = lm(Y~condition, data = data_file_sub)
            out = summary(res)$coefficients[2,c(1,4)]
            names(out) = c(est_head, "pvalue")            
        } else if (test == "m") {
            res = lme4::lmer(Y~condition + (1 | lib_id), data = data_file_sub)
            out = summary(res)$coefficients[2,c(1,3)]
            out[2] = 2*pt(abs(out[2]), df=nrow(data_file_sub)-2,lower.tail = FALSE)
            names(out) = c(est_head, "pvalue")
        }}
    ,
        error = {function(e) {warning(
            sprintf("Error: NAs given"))
        }}
    )
    
    out$qvalue = p.adjust(out$pvalue, n = "BH")

    return(out)
}

diff_length = function(data_file, test) {
    #genes = names(table(data_file$id))
    
    data_file_bygene = split(data_file, data_file[,2])
    
    out = lapply(data_file_bygene, function(d) {diff_length_single(d, test)})
    out = do.call(rbind, out)
    return(out)
}

outres = diff_length(data_file, test)

write.table(outres, file = ofile,sep = delim, quote = F, row.names = T, col.names = T)
