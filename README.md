## Introduction

This code was provided to me by Sarah Pyfrom.



## Installation

This package was created by running the following in an R terminal:

```R
install.packages("devtools")
install.packages("seqinr")
install.packages("logr")

install.packages("BiocManager")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")

devtools::create("/home/harrisonized/github/R/escapegenecalculator")
devtools::document("/home/harrisonized/github/R/escapegenecalculator")
```



## Downloading Data

Check the README files for what to download for each folder.



## Use

Run the following from command line:

```bash
cd escapegenecalculator
Rscript R/calculate_sequence_lengths.R
```

