## Note: this script takes about 20 seconds to run
## 1. Import the gzipped gtf file using rtracklayer
## 2. convert it to a data frame
## 3. throw away rows that have an NA in the gene_name column
## 4. throw out rows that have genes mapped to autosomes
## 5. pull out the chromosome, start, stop, width, type, gene id, and gene name columns.
## 6. Loop through and pull out the width and gene names that correspond to type == 'gene'.
## 7. Do for both genomes.


library(logr)
library(GenomicFeatures)
wd = dirname(this.path::here())


# Input parameters
# Note: in the future, make this loop through the directory or be able to pass a list
in_dir = file.path(wd, "data/ref/gtf")
in_filename = 'Mus_musculus_casteij.CAST_EiJ_v1.103.chr_sorted.gtf.gz'
# in_filename = 'Mus_musculus.GRCm38.102.chr_sorted_noYMT.gtf.gz'
out_dir = file.path(wd, "data/ref")
species = strsplit(in_filename, split = ".", fixed=TRUE)[[1]][1]
out_filename = paste("exon_lengths-", species, ".csv", sep='')


# Start Log
log <- log_open(paste("calculate_exon_lengths ", Sys.time(), '.log', sep='')) # Open log
log_print(paste('input file: ', file.path(in_dir, in_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

# Import gtf data
log_print("Importing gtf data...")
# gtf <- rtracklayer::import('data/gtf/Mus_musculus.GRCm38.87.gtf')
# gtf <- rtracklayer::import('data/gtf/Mus_musculuseij.CAST_EiJ_v1.108.gtf')
gtf_data <- rtracklayer::import(
    file.path(in_dir, in_filename)
)
gtf_df <- as.data.frame(gtf_data, stringsAsFactors = FALSE)
gtf_df_x_only <- gtf_df[is.na(gtf_df$gene_name) == FALSE & gtf_df$seqnames == 'X',]
# gtf_subset <- gtf_df_x_only[,c(1,2,3,4,7,10,23)]
# gtf_subset <- gtf_df_x_only[,c(1,2,3,4,7,10,15)]  # this is bad
gtf_subset <- gtf_df_x_only[, c("seqnames", "start", "end",  "width", "type", "gene_id", "gene_name")]


# First, import the GTF-file that you have also used as input for htseq-count
start_time = Sys.time()
log_print(paste("Making TxDb from GTF...", start_time))
txdb <- makeTxDbFromGFF(
    file.path(in_dir, in_filename),
    format="gtf"
)
end_time = Sys.time()
log_print(paste('finished.', end_time))
log_print(paste('time of operation... ', end_time - start_time))


# ----------------------------------------------------------------------
# Processing

# then collect the exons per gene id
exons.list.per.gene.cast <- exonsBy(txdb,by="gene")

# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes.cast <- as.data.frame(sum(width(reduce(exons.list.per.gene.cast))))

# try to match exonic lengths to genes on the X using the gtf_subset df and exonic.gene.sizes df.
# final_data <- gtf_subset[,c(6,7)]
final_data <- gtf_subset[,c("gene_id", "gene_name")]  # Pulling out just gene_ids and gene_names from the gtf file
final_data <- unique(final_data)  # final_dataping down this new df to only unique instances

index_in_exonic.gene.sizes.cast <- match(final_data[,c("gene_id")], rownames(exonic.gene.sizes.cast))  # find the indexes of the entries in exonic.gene.sizes that matches the genes in the final_data dataframe.

final_data[,3] <- exonic.gene.sizes.cast[index_in_exonic.gene.sizes.cast,] # Append those indexes onto the final_data df and rename cols. 
colnames(final_data) <- c('gene_id', 'gene_name', 'exon_length')


# ----------------------------------------------------------------------
# Save Output

log_print(paste('writing file...', Sys.time()))
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

log_print(paste('End', Sys.time()))
log_close()