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
        } else if (test == "s") {
          library(splines)
          if (is.null(params)) {
            model0 = length ~ (1 | lib_id)
            model1 = length ~ (1 | lib_id) + bs(condition, degree = 2)
          } else {
            model0 = as.formula(sprintf("length ~ %s + (1 | lib_id)", params))
            model1 = as.formula(sprintf("length ~ %s + (1 | lib_id) + bs(condition, degree = 2)", params))
          }
          
          suppressMessages({
            f0 = lme4::lmer(model0, data = data_file_sub, REML = F)
            f1 = lme4::lmer(model1, data = data_file_sub, REML = F)
          })
          out = data.frame(p = anova(f0,f1)$`Pr(>Chisq)`[2])
          est_head = NULL
        }
      }
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
    suppressWarnings({out = data.frame(do.call(rbind, out))})
    rownames(out) = names(data_file_byname)
    out$qvalue = p.adjust(out$pvalue,method = "BH")
    out = out[!is.na(out$pvalue), ]
    desc = lapply(rownames(out),
                  function(d) {calc_descriptives(data_file_byname[[d]])})
    desc = data.frame(do.call(rbind, desc))
    if (!is.null(baseline)) {
      colnames(desc) = c(paste("n",b,sep = "."),"n.alt",
                         paste("mean_length",b,sep = "."),"mean_length.alt")
    } else {
      conds = names(table(data_file$condition))
      colnames(desc) = c(paste("n.", conds, sep = ""), paste("mean_length.", conds, sep = ""))
    }
    out = cbind(out, desc)
    return(out)
}
