## Introduction

This code was provided to me by Sarah Pyfrom. It takes tsv files from `data/raw_read_counts` as input and outputs a variety of csv files.



## Installation

Run the following in your R console. My machine is on R 4.2.2.

```R
install.packages("devtools")
install.packages("seqinr")
install.packages("logr")
install.packages("data.table")
install.packages("dplyr")

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

# generate intermediate files for use in step1
Rscript R/calculate_exon_lengths.R
Rscript R/calculate_chromosome_lengths.R  # this takes 5 minutes

# Go through the steps
Rscript R/step0-merge_read_counts.R
Rscript R/step1-calculate_RPKMs_SRPMs_from_RPMs.R
Rscript R/step2-aggregate_raw_reads.R
Rscript R/step3-CI_calculation.R
Rscript R/step4-match_CIs_RPKM_SRPM.R
```



## Style Guide

1. Use snake_case, not camelCase.
2. Make sure all inputs/outputs are csv, not tsv.