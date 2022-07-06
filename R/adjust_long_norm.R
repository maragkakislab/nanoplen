
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

run_adjust_norm = function(data_file, metadata) {
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
    
    return(data_file)
}


