#' Run NanopLen differential length on a file of lengths per transcript
#'
#' NanopLen tests for significant differences in lengths, e.g. transcript length
#' or Poly(A) tail length, between two categories. 
#' 
#' @section Expected inputs:
#' \describe{
#'   \item {data file}{A file of lengths with these columns in order: library id, gene/transcript id, length}
#'   \item {metadata}{Metadata with these columns in order: library id, condition, and any optional additional columns}
#'   \item {test}{Type of statistical test: 't' for t-test, 'm' for linear mixed model, 'w' for Wilcoxon}
#'   \item {baseline}{String that indicates the baseline variable for the 'condition' variable in the metadata}
#'   \item {logscale}{ indicates if length should be on log2scale}
#'   \item {params}{Extra parameters to use in the model, done in model notation (e.g. to add age and sex and their interaction to the model, the input is 'age+sex+age*sex') }
#' }
#' 
#' @section Expected output:
#' A results file with the following columns:
#' \describe{
#'   \item {id}{corresponds to gene/transcript as given in data file}
#'   \item {meandiff}{the mean difference in length between conditions, shown if logscale = F}
#'   \item {log2fc}{log2 fold change in length between conditions, shown if logscale = T}
#'   \item {pvalue}{p-values for significance of differential length test}
#'   \item {qvalue}{Adjusted p-values using the Bonferroni-Hochberg method}
#'   \item {n.baseline}{Number of observations for baseline condition. "baseline" will be replaced with user-defined baseline string}
#'   \item {n.alt} {Number of observations for alternate condition}
#'   \item {mean_length.baseline}{Mean length of reads in baseline condition. "baseline" will be replaced with user-defined baseline string}
#'   \item {mean_length.baseline}{Mean length of reads in alternate condition.}
#' }
#' 
#' @section Tests:
#' \describe{
#'   \item {t}{The most basic t-test between conditions. Supports additional parameters which makes this functionally equivalent to
#' simple linear regression}
#'   \item {m}{Linear Mixed Model, testing for difference between conditions with a random effect on the libraries. Supports additional parameters.}
#'   \item {w}{Wilcoxon Rank-Sum test, testing for difference between conditions. Does not support extra parameters.}
#' }
#' 
#' @export
#' @include diff_length_funcs.R

nanoplen <- function(data_file,
                    metadata,
                    test = "t",
                    baseline = "Control",
                    logscale = FALSE,
                    params = NULL,
                    norm = FALSE
) {
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
    
    # Remove rows with reported length 0. They should not be there anyway.
    data_file = data_file[data_file$length > 0,]
    
    # If normalizing, check if norm_group column is in metadata
    if (norm) {
        if (!("norm_group" %in% colnames(metadata))) {
            stop("Need 'norm_group' column in metadata!")
        }
        data_file = run_adjust_norm(data_file, metadata)
    }
    
    # Add condition column from metadata
    data_file = merge(data_file, metadata, by="lib_id")[,1:4]
    
    if (test == "w") {
        levels = levels(metadata$condition)
        logscale = FALSE
        if (length(levels)>2) {warning("More than two levels detected, Wilcox is only for two-level comparison!\n")}
        if (!is.null(params)) {warning("Extra parameters not supported with Wilcoxon test!\n")}
    }
    
    # Relevel data_file$condition to use baseline string
    data_file = within(data_file, condition <- relevel(factor(condition), ref = baseline))
    
    outres = diff_length(data_file, test, params, logscale, baseline)
    outres = cbind(rownames(outres),outres)
    colnames(outres)[1] = "name"
    
    return(outres)
}
