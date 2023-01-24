## Introduction

This code was provided to me by Sarah Pyfrom. Original author was Zachary Beetham. It takes tsv files from `data/read_counts` as input and outputs a variety of csv files. Implementation is based on this paper by Joel Berletch et al.:

```
Berletch JB, Ma W, Yang F, Shendure J, Noble WS, Disteche CM, et al. (2015) Escape from X Inactivation Varies in Mouse Tissues. PLoS Genet 11(3): e1005079. doi:10.1371/journal.pgen.1005079
```



## Installation

Run the following in your R console. My machine is on R 4.2.2.

```R
install.packages("devtools")
install.packages("seqinr")
install.packages("logr")
install.packages("data.table")
install.packages("dplyr")
install.packages("this.path")  # See: https://github.com/ArcadeAntics/this.path

install.packages("BiocManager")
BiocManager::install("rtracklayer")  # this takes ~15 minutes
BiocManager::install("GenomicFeatures")

devtools::create("~/github/R/escapegenecalculator")
devtools::document("~/github/R/escapegenecalculator")

# Warning! Do not install this as a library, it will run through all scripts
# install.packages("~/github/R/escapegenecalculator", repos=NULL, type='source')  
```



## Downloading Data

Check the README files for what to download for each folder.



## Use

Run the following from command line:

```bash
cd escapegenecalculator

# generate metadata files for use downstream
Rscript R/calculate_exon_lengths.R  # ~1 minute, required
Rscript R/extract_chromosome_lengths.R  # ~5 minutes, optional

# Go through the steps
Rscript R/merge_read_counts.R
Rscript R/step1-calculate_confidence_intervals.R
Rscript R/step2-calculate_rpkm_srpm.R
Rscript R/step3-filter_rpkm_sprm_again  # will be deleting this
```



## Style Guide

1. Use snake_case, not camelCase.
2. Make sure all outputs are csv, **NOT** tsv.