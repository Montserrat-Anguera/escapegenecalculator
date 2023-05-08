## Introduction

The purpose of this R repository is to calculate escape genes. Implementation is based on this paper by Joel Berletch et al.:

```
Berletch JB, Ma W, Yang F, Shendure J, Noble WS, Disteche CM, et al. (2015)
Escape from X Inactivation Varies in Mouse Tissues. PLoS Genet 11(3): e1005079.
doi:10.1371/journal.pgen.1005079
```

This code has existed in Montserrat Anguera's lab for a while, passing through multiple hands before it finally made its way to me (Harrison Wang). Xiang Yu wrote the original code some time in 2020. It was then passed to Zachary Beethem, who maintained and modified it until he left in early 2021, after which it was passed to Sarah Pyfrom. I received these files from Sarah.


While working on this project, one of my early questions is why our lab consistently comes up with many more escape genes (~25%) than what the literature suggests (<10%). After much testing, finally going back to Berletch's source data, located on GEO here: [GSE59777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59777), we finally discovered three major errors:

1. Zack's version systematically overestimated the RPKMs (off by ~2 orders of magnitude) because it was systematically underestimating the `total_num_reads`, as it misinterpreted that quantity to be the sum of SNP-specific reads. In fact, the correct `total_num_reads` must be calculated on all mapped reads, not just SNP-specific reads, so it cannot be computed from the input into this pipeline. Rather, `total_num_reads` must be supplied a priori. 
2. Zack's version averages mice together, but in Disteche's chapter in the X Chromosome Inactivation Methods book published in 2018 (found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6269188/)), they say this: "Biological replicates of RNA-seq experiments should be analyzed separately." 
3. Finally, like the `total_num_reads`, the RPKMs must also be derived a priori, because they are calculated on all mapped reads, not just the SNP-specific reads. In the absence of that information, I impute the `total_num_reads` by estimating that the total SNP-specific reads is about 0.121 of the `total_num_reads`. Using this as a estimate brings the computed/estimated RPKM down to reasonable values.

With these changes, we now get a much shorter list of escape genes, more consistent with what the literature suggests.


## Installation

Run the following in your R console. My machine is on R 4.3.0, but I believe any version should work fine, as I'm only using Base R and nothing fancy.

```R
install.packages("devtools") # warning: takes a long time!
install.packages("reshape")
install.packages("seqinr")
install.packages("logr")
install.packages("data.table")
install.packages("this.path")  # See: https://github.com/ArcadeAntics/this.path
install.packages("optparse")
install.packages("tidyr")
# install.packages("dplyr")  # eliminated this as a dependency

install.packages("BiocManager")
BiocManager::install("rtracklayer")  # this takes ~15 minutes
BiocManager::install("GenomicFeatures")

# documentation only:
devtools::create("~/path/to/escapegenecalculator")
devtools::document("~/path/to/escapegenecalculator")

# Warning! Do not install this as a library, it will run through all scripts
# install.packages("~/github/R/escapegenecalculator", repos=NULL, type='source')  
```



## Downloading Data

Check the README files for what to download for each folder.



## Use

All the combinations are listed here for your convenience.

Preprocessing

```bash
cd path/to/escapegenecalculator

# Optional, since the output of these have already been generated for you
Rscript R/calculate_exon_lengths.R  # ~1 minute
Rscript R/extract_chromosome_lengths.R  # ~5 minutes
```

Old pipeline (output in `output-1`):

```bash
cd path/to/escapegenecalculator

# data/berletch-spleen
Rscript R/generate_config.R -i data/berletch-spleen -o output-1
Rscript R/escapegenecalculator_old.R -i data/berletch-spleen

# data/sierra-at2
Rscript R/generate_config.R -i data/sierra-at2 -o output-1
# make sure to update the mouse_gender
Rscript R/escapegenecalculator_old.R -i data/sierra-at2
```

New pipeline (output in `output-2`):

```bash
cd path/to/escapegenecalculator

# data/berletch-spleen
Rscript R/generate_config.R -i data/berletch-spleen -o output-2 # default
Rscript R/escapegenecalculator.R -i data/berletch-spleen

# data/sierra-at2
Rscript R/generate_config.R -i data/sierra-at2 -o output-2
Rscript R/escapegenecalculator.R -i data/sierra-at2
```

## Contributing

Reach out to harrison.c.wang@gmail.com and we can set up a time to discuss. You could also make your own github branch or file a github issue.