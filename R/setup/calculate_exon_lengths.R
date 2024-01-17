## This script takes about 20 seconds to run
## 1. Import the gzipped gtf file using rtracklayer
## 2. convert it to a data frame
## 3. throw away rows that have an NA in the gene_name column
## 4. throw out rows that have genes mapped to autosomes
## 5. pull out the chromosome, start, stop, width, type, gene id, and gene name columns.
## 6. Loop through and pull out the width and gene names that correspond to type == 'gene'.
## 7. Do for both genomes.

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('GenomicFeatures')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-file"),
                default="ref/gtf/Mus_musculus_casteij.CAST_EiJ_v1.103.chr_sorted.gtf.gz",
                # default = "ref/gtf/Mus_musculus.GRCm38.102.chr_sorted_noYMT.gtf.gz",
                metavar="ref/gtf/Mus_musculus_casteij.CAST_EiJ_v1.103.chr_sorted.gtf.gz", type="character",
                help="set the input file"),

    make_option(c("-o", "--output-dir"), default="ref",
                metavar="ref", type="character",
                help="data will be output in this folder inside the input-dir"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

input_file = file.path(wd, opt[['input-file']])
species = strsplit(basename(input_file), split = ".", fixed=TRUE)[[1]][1]
out_dir = file.path(wd, opt[['output-dir']])
out_filename = paste("exon_lengths-", species, ".csv", sep='')


# Start Log
start_time = Sys.time()
log <- log_open(paste0("calculate_exon_lengths-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('input file: ', input_file))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read and Preprocess Data

# Import gtf data
log_print("Importing gtf data...")
gtf_data <- rtracklayer::import(file.path(wd, opt[['input-file']]))
# gtf_data <- rtracklayer::import('data/gtf/Mus_musculus.GRCm38.87.gtf')
# gtf_data <- rtracklayer::import('data/gtf/Mus_musculuseij.CAST_EiJ_v1.108.gtf')

# Preprocessing
gtf_df <- as.data.frame(gtf_data, stringsAsFactors = FALSE)
gtf_df_x_only <- gtf_df[is.na(gtf_df$gene_name) == FALSE & gtf_df$seqnames == 'X',]
gtf_subset <- gtf_df_x_only[, c("seqnames", "start", "end",  "width", "type", "gene_id", "gene_name")]

# Deprecated:
# gtf_subset <- gtf_df_x_only[,c(1,2,3,4,7,10,23)]
# gtf_subset <- gtf_df_x_only[,c(1,2,3,4,7,10,15)]


# First, import the GTF-file that you have also used as input for htseq-count
process_start_time = Sys.time()
log_print(paste("Making TxDb from GTF...", process_start_time))

txdb <- makeTxDbFromGFF(
    file.path(in_dir, in_filename),
    format="gtf"
)

process_end_time = Sys.time()
log_print(paste('finished.', process_end_time))
log_print(paste('time of operation... ', process_end_time - process_start_time))


# Collect the exons per gene id
exons.list.per.gene.cast <- exonsBy(txdb,by="gene")

# For each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes.cast <- as.data.frame(sum(width(reduce(exons.list.per.gene.cast))))

# Try to match exonic lengths to genes on the X using the gtf_subset df and exonic.gene.sizes df.
# final_data <- gtf_subset[,c(6,7)]
final_data <- gtf_subset[,c("gene_id", "gene_name")]  # Pulling out just gene_ids and gene_names from the gtf file
final_data <- unique(final_data)  # final_dataping down this new df to only unique instances

# find the indexes of the entries in exonic.gene.sizes that matches the genes in the final_data dataframe.
index_in_exonic.gene.sizes.cast <- match(final_data[,c("gene_id")], rownames(exonic.gene.sizes.cast)) 

# Append those indexes onto the final_data df and rename cols
final_data[,3] <- exonic.gene.sizes.cast[index_in_exonic.gene.sizes.cast,]
colnames(final_data) <- c('gene_id', 'gene_name', 'exon_length')


# save
if (!troubleshooting) {

    log_print(paste('Writing file...', Sys.time()))

    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }

    write.table(
        final_data,
        file = file.path(out_dir, out_filename),
        col.names = TRUE,
        row.names = FALSE,
        sep = ','
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
