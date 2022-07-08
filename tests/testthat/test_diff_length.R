context("Test diff_length functions")

#Enter R in package root, run `test_local()`

# test_that('label' , {
#   [execute code]
#   expect_*(...)
# } )

testdata = read.table("../test_data/plen_test_data.tab", sep="\t",header=T)
metadata = read.table("../test_data/plen_test_metadata.tab", sep="\t",header=T)
colnames(testdata) = c("lib_id","name","length")
colnames(metadata)[1:2] = c("lib_id", "condition")
data_file = merge(testdata, metadata, by="lib_id")[,1:4]

testdata_bad = read.table("../test_data/plen_test_data_bad.tab", sep="\t",header=T)
colnames(testdata_bad) = c("lib_id","name","length")
data_file_bad = merge(testdata_bad, metadata, by="lib_id")[,1:4]

test_that('Descriptives', {
    df = data.frame(length = 1:10, condition = rep(c("cond","ctrl"),5))
    out = calc_descriptives(df)
    expect_equal(out, c(5,5,5,6))
})

test_that('test_single_t', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "t")
    expect_equal(names(out), c("log2FC","pvalue"))
})

test_that('test_single_t_nolog', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "t", logscale = FALSE)
    expect_equal(names(out), c("meandiff","pvalue"))
})

test_that('test_single_t_extravars', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "t",params = "zee")
    expect_equal(names(out), c("log2FC","pvalue"))
})

test_that('test_single_mix', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "m")
    expect_equal(names(out), c("log2FC","pvalue"))
})

test_that('test_single_mix_nolog', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "m", logscale = FALSE)
    expect_equal(names(out), c("meandiff","pvalue"))
})

test_that('test_single_w', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "w")
    expect_equal(names(out), c("Wilcox_stat","log2FC","pvalue"))
})

test_that('test_single_w_nolog', {
    df = data_file[data_file$name == "gene_1",]
    out = diff_length_single(df, "w", logscale = FALSE)
    expect_equal(names(out), c("Wilcox_stat","log2FC","pvalue"))
})

test_that('test_single_bad', {
    df = data_file_bad[data_file_bad$name == "gene_3",]
    expect_error(out = diff_length_single(df, "t"))
}) 

test_that('all_test_t', {
    out = diff_length(data_file, "t", NULL, TRUE, "control")
    expect_equal(nrow(out), 2)
    expect_equal(colnames(out), c("log2FC", "pvalue","qvalue","n.control","n.alt","mean_length.control","mean_length.alt"))
})

test_that('all_test_mix', {
    out = diff_length(data_file, "m", NULL, TRUE, "control")
    expect_equal(nrow(out), 2)
    expect_equal(colnames(out), c("log2FC", "pvalue","qvalue","n.control","n.alt","mean_length.control","mean_length.alt"))
})

test_that('all_test_w', {
    out = diff_length(data_file, "w", NULL, TRUE, "control")
    expect_equal(nrow(out), 2)
    expect_equal(colnames(out), c("Wilcox_stat","log2FC", "pvalue","qvalue","n.control","n.alt","mean_length.control","mean_length.alt"))
})

test_that('all_test_bad', {
    expect_warning({out = diff_length(data_file_bad, "t", NULL, TRUE, "control")})
    expect_equal(nrow(out), 3)
    expect_true(is.na(out[3,1]))
})


test_that('test_single_t_value', {
    df = data_file[data_file$name == "gene_2",]
    out = diff_length_single(df, "t", logscale = FALSE)
    expect_true(round(out["meandiff"],0)==-43 )
})

test_that('test_single_t_value_reversebaseline', {
    data_file_2 = within(data_file, condition <- relevel(condition, ref = "treated"))
    df = data_file_2[data_file_2$name == "gene_2",]
    out = diff_length_single(df, "t", logscale = FALSE, )
    expect_true(round(out["meandiff"],0)==43 )
})