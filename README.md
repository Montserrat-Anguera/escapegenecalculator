## Introduction

This code was provided to me by Sarah Pyfrom.



## Installation

Run the following in your R console. My machine is on R 4.2.2.

```R
install.packages("devtools")
install.packages("seqinr")
install.packages("logr")

install.packages("BiocManager")
BiocManager::install("rtracklayer")  # this takes ~15 minutes
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

# generate intermediate files for use in step1
Rscript R/calculate_exon_lengths.R
Rscript R/calculate_chromosome_lengths.R  # this takes 5 minutes

# Go through the steps
Rscript R/step1-calculate_RPKMs_SRPMs_from_RPMs.R
Rscript R/step2-convert_rawreads_to_filteredreads.R
Rscript R/step3-CI_calculation.R
Rscript R/step4-match_CIs_RPKM_SRPM.R
```
