diff_length_single = function(data_file_sub, test, params = NULL, logscale = TRUE, contig_name = "") {
    out = c(NA,NA)
    if (is.null(params)) {
        model = "length~condition"
    } else {
        model = paste(c("length~condition",params),sep="+")
    }
    if (logscale) {
        data_file_sub$length = log2(data_file_sub$length)
        est_head = "log2FC"
    } else {
        est_head = "meandiff"
    }
    tryCatch({   
        if (test == "t") {
            res = lm(as.formula(model), data = data_file_sub)
            out = summary(res)$coefficients[2,c(1,4)]
        } else if (test == "m") {
            model = paste(model, "(1 | lib_id)", sep = "+")
            res = lme4::lmer(as.formula(model), data = data_file_sub)
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
            sprintf("Error in %s: NAs given\n", contig_name))
            has_warning <<- TRUE 
        }}
    )
    names(out) = c(est_head, "pvalue")
    
    return(out)
}

calc_descriptives = function(df) {
    out = c(NA,NA,NA,NA)
    
    tryCatch({  
        out = c(aggregate(df$length, list(df$condition), FUN = length)[,2],
                aggregate(df$length, list(df$condition), FUN = mean)[,2])
    }, error = function(e) { }
    )
    return(out)
}

diff_length = function(data_file, test, params, logscale, b = baseline) {
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
