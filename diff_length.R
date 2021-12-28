#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

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
                help="Normalize data by replicate. Requires replicate column in metadata where replicates have the same value with each other [default %default]"),
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

# Can omit this step if previous output names columns consistently
colnames(data_file) = c("lib_id","name","length")
colnames(metadata)[1:2] = c("lib_id", "condition")

# checking if model parameters are in metadata
if (!is.null(params)) {
    vars = unique(unlist(strsplit(strsplit(params,"\\+")[[1]], "\\*")))
    vars_in_meta = vars %in% colnames(metadata)
    if (!all(vars_in_meta)) {
        stop(sprintf("Model parameters not in metadata: %s",paste(vars[!vars_in_meta], collapse = " ")))
    }
}

# If normalizing, check if replicate column is in metadata
if (norm) {
    if (!("replicate" %in% colnames(metadata))) {
        stop("Need 'replicate' column in metadata!")
    }
}

# Remove rows with reported length 0. They should not be there anyway.
data_file = data_file[data_file$length > 0,]

# Normalize data by shifting the longer reads to match 1-1 correlation
# Returns both columns with the shorter normalized
adjust_long_norm = function(x,y) {
    t1 = prcomp(cbind(x,y))
    # Slope
    b11 = t1$rotation[2,1]/t1$rotation[1,1]
    # Intercept
    b10 = mean(y) - b11*mean(x)
    
    c = b10/(1-b11)
    
    d = cbind(x,y)
    #d new
    dn = d
    
    if (b11>1) {
        dn[dn[,2]>c,1] = dn[dn[,2]>c,1]+dn[dn[,2]>c,2]+b10/b11-1/b11*dn[dn[,2]>c,2]
    } else {
        dn[dn[,1]>c,2] = dn[dn[,1]>c,2]+dn[dn[,1]>c,1]-b10-b11*dn[dn[,1]>c,1]
        
    }
    
    return(dn)
    
}

adjust_long_norm_many = function(data_sub) {
    n = ncol(data_sub)
    if (n == 1) {
        # No need to normalize
        warning(sprintf("%s has no replicates!",))
        return(data_sub)
    }
    
    # Make an upper triangular comparison
    pcomps = matrix(0, n,n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                pcomps[i,j] = 1 #to be doubled during the mirroring step
            } else if (i < j) {
                t1 = stats::prcomp(cbind(data_sub[,i],data_sub[,j]))
                pcomps[i,j] = t1$rotation[2,1]/t1$rotation[1,1]
            } else {
                pcomps[i,j] = 1/pcomps[j,i]
            }
        }
    }
    
    # Find out the library with the steepest slope
    top = which(sapply(1:n, function(i) {all(pcomps[i,]>=1)}))[1]
    
    data_sub_norm = data_sub
    for (i in 1:n) {
        if (i != top) {
            data_sub_norm[,c(i,top)] = adjust_long_norm(data_sub[,i],data_sub[,top])
        }
    }
    
    colnames(data_sub_norm) = colnames(data_sub)
    return(data_sub_norm)    
}

if (norm) {
    rep_ind = split(1:nrow(metadata), metadata$replicate)
    
    #Convert from long to wide, aggregating mean lengths
    data_mean_length = aggregate(length ~ lib_id+name, data = data_file, FUN = mean) 
    data_length_wide = reshape(data_mean_length, idvar = "name", timevar = "lib_id", direction = "wide")
    #Remove transcripts where there are NAs
    data_length_wide = data_length_wide[!apply(data_length_wide, 1, 
                                              function(t) {any(is.na(t))}),]
    
    rownames(data_length_wide) = data_length_wide$name
    data_length_wide = data_length_wide[,-1]
    colnames(data_length_wide) = sapply(colnames(data_length_wide), function(x) {substr(x, 8,1E5)})
    data_length_wide = data_length_wide[,metadata$lib_id]
    
    
    data_length_wide = log2(data_length_wide)
    length_adjust = data_length_wide
    for (inds in rep_ind) {
        length_adjust[,inds] = adjust_long_norm_many(data_length_wide[,inds])
    }
    
    # Find difference matrix
    length_adjust = length_adjust - data_length_wide
    
    # Adjust long format
    length_adjust = 2^(length_adjust)
    length_adjust$name = rownames(length_adjust)
    length_adjust_long = reshape(length_adjust, direction = "long", idvar = "name",
                                 varying = colnames(length_adjust)[-ncol(length_adjust)],
                                 times = colnames(length_adjust)[-ncol(length_adjust)],
                                 v.names = "length_adj",
                                 timevar = "lib_id")
    
    data_file = dplyr::inner_join(data_file, length_adjust_long, by = c("name", "lib_id"))
    data_file$length = data_file$length * data_file$length_adj
    data_file = data_file[,c("lib_id","name","length")]

}

# Add condition column from metadata
data_file = merge(data_file, metadata, by="lib_id")[,1:4]

if (test == "w") {
    levels = levels(metadata$condition)
    logscale = FALSE
    if (length(levels)>2) {warning("More than two levels detected, Wilcox is only for two-level comparison!")}
    if (!is.null(params)) {warning("Extra parameters not supported with Wilcoxon test!")}
}

# Relevel data_file$condition to use baseline string
data_file = within(data_file, condition <- relevel(condition, ref = baseline))

has_warning <<- FALSE
wfile = ifelse(ofile == "stdout", "warnings.txt", sprintf("%s_warnings.txt", ofile))
ww <- file(wfile, open = "wt")
sink(ww, type = "message")

diff_length_single = function(data_file_sub, test, params, logscale = TRUE, contig_name) {
    out = c(NA,NA)
    model = paste(c("length~condition",params),sep="+")
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
        } else if (test == "w") {
            levels = levels(data_file_sub$condition)
            x = data_file_sub$length[data_file_sub$condition == levels[1]]
            # Supposed to only have two levels, but just in case
            y = data_file_sub$length[data_file_sub$condition != levels[1]]
            # Code from wilcox.test.default
            r <- rank(c(x, y))
            n.x <- as.double(length(x))
            n.y <- as.double(length(y))
            STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
            NTIES <- table(r)
            z <- STATISTIC - n.x * n.y/2
            SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) - 
                                                sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 
                                                                                         1))))
            z = z/SIGMA
            r_pval = 2 * min(stats::pnorm(z),
                             stats::pnorm(z, lower.tail = FALSE))  
            l2f = log2(mean(y)/mean(x))
            out = c(STATISTIC, l2f, r_pval)
            est_head = c("Wilcox_stat","log2FC")
        }}
        ,
        error = {function(e) {warning(
            #FIXME: Was supposed to also show which gene/transcript but cannot extract with current algorithm
            sprintf("Error in %s: NAs given", contig_name))
            has_warning <<- TRUE 
        }}
    )
    names(out) = c(est_head, "pvalue")
    
    return(out)
}

calc_descriptives = function(d) {
    out = c(NA,NA,NA,NA)
    
    tryCatch({  
        out = c(aggregate(d$length, list(d$condition), FUN = length)[,2],
                aggregate(d$length, list(d$condition), FUN = mean)[,2])
    }, error = function(e) { }
    )
    return(out)
}

diff_length = function(data_file, test, params, b = baseline) {
    data_file_byname = split(data_file, data_file$name)
    
    # Loops over all subsets split by name
    out = lapply(names(data_file_byname), 
                 function(d) {diff_length_single(
                     data_file_byname[[d]], test, params, logscale, d)})
    out = data.frame(do.call(rbind, out))
    out$qvalue = p.adjust(out$pvalue,method = "BH")
    desc = lapply(names(data_file_byname),
                  function(d) {calc_descriptives(data_file_byname[[d]])})
    desc = data.frame(do.call(rbind, desc))
    colnames(desc) = c(paste("n",b,sep = "."),"n.alt",
                       paste("mean_length",b,sep = "."),"mean_length.alt")
    out = cbind(out, desc)
    rownames(out) = names(data_file_byname)
    return(out)
}

outres = diff_length(data_file, test, params)
outres = cbind(rownames(outres),outres)
colnames(outres)[1] = "name"

if (ofile == "stdout") {
    write.table(outres,file=stdout(),sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
} else {
    write.table(outres, file = ofile,sep = delim, quote = F, row.names = FALSE, col.names = TRUE)
}

sink(type="message")
close(ww)
if (has_warning) {
    message(sprintf("Some tests errored, logged in %s", wfile))
}
