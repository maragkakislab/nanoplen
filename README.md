# nanoplen: Differential length analysis
Definition: Nanoplen tests for differences in transcript length between diffrerent conditions in nanopore sequencing data. 

# Nanoplen differential length model:
Nanoplen has the option to use three models: t-test, linear mixed model, or Wilcoxon test. Regardless of test used, the output is the pertaining statistic and p-value for each contig in the data. T-test and linear mixed model can include extra covariates but Wilcoxon cannot.
t-test:
Nanoplen models a simple linear regression for the effect of the condition on average length. If no extra covariates are supplied, this is mathematically equivalent to a t-test for difference in means.
Linear Mixed Model:
The linear mixed model uses the library ID as a random effect to control for library-specific batch effects. Otherwise, this uses the same model as in the t-test.
Wilcoxon:
The Wilcoxon test is the nonparametric equivalent of the t-test for difference in means. This method is restricted to a two-level condition variable and no extra covariates.

we used LMM over t-test because it can account for sample-specific batch effects and is just more powerful in general. LMM was also used over wilcoxon for similar reasons, but also because wilcoxon cannot account for extra parameters, but we just put it there in case the length distribution is all over the place for some reason


# Nanoplen normalization:
We assume that discrepancies in length in replicates is far more likely to be shortening of longer transcripts, while shorter transcripts are unaffected. To normalize two replcates to each other, we first use a simple linear regression model on the average length of each transcript in the two replicates. We use the slope to determine which replicate has larger average transcript lengths. We then increase the length of all transcripts above the intersection between the modeled regression line and the Y=X of the sample with smaller average lengths to match the Y=X line. For normalization with more than two replicates, all replicates are adjust to match the replicate with the longest transcripts.

To use this option, the metadata must require a `norm_group` column, where replicates have identical values.

## Input

Nanoplen accepts two input files. The first file is TAB separated and contains the length data and consists of 3 columns, the sample name, the identifier and the length. The second file is TAB separated and contains the metadata that describe the samples and the conditions in the experiment. The identifier is used to aggregate lengths together, usually a gene or transcript id.

### Example inputs

Length data file called lengths.txt

```
lib_id  id  len
library_1 gene_1  1500
library_1 gene_1  1400
library_1 gene_1  1450
library_1 gene_2  500
library_1 gene_2  600
library_2 gene_1  1600
library_2 gene_2  650
library_3 gene_1  1200
library_3 gene_1  1100
library_3 gene_1  1300
library_3 gene_2  400
library_4 gene_1  1000
library_4 gene_2  800
library_4 gene_2  500
library_5 gene_1  1350
library_5 gene_1  1300
library_5 gene_1  1200
library_5 gene_2  200
library_5 gene_2  800
```

Metadata file, called metadata.txt

```
lib_id  condition
library_1 control
library_2 control
library_3 treat
library_4 treat
library_5 treat
library_6 treat
```


### Example command
```
diff_length.R -d lengths.txt -m metadata.txt -t w --logscale -o path/to/out
```

### Example output
```
id  baseMean  log2fc  pvalue
gene_1  1500  -0.5  0.0001
gene_2  500 -0.1  0.1
```
