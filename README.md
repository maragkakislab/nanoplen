# nanoplen: Differential length analysis

## Input

Nanoplen accepts two input files. The first file is TAB separated and contains the metadata that describe the samples
and the conditions in the experiment, similar to DESeq2. The second file is also TAB separated and contains the actual
length data and consists of 3 columns, the sample name, the identifier and the length. The identifier is used to
aggregate lengths together. It is usually a gene or transcript id but it can be anything else that the user wishes.

### Example input

Metadata file

```
lib_id  condition
library_1 control
library_2 control
library_3 treated
library_4 treated
library_3 treated
library_5 treated
```

Length data file
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

### Example output
```
id  baseMean  log2fc  pvalue
gene_1  1500  -0.5  0.0001
gene_2  500 -0.1  0.1
```
