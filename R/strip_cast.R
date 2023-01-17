## 1. Import the gzipped gtf file using rtracklayer
## 2. convert it to a data frame
## 3. throw away rows that have an NA in the gene_name column
## 4. throw out rows that have genes mapped to autosomes
## 5. pull out the chromosome, start, stop, width, type, gene id, and gene name columns.
## 6. Loop through and pull out the width and gene names that correspond to type == 'gene'.
## 7. Do for both genomes. 


library(logr)
library(GenomicFeatures)


# Input parameters
gtf_dir = file.path(getwd( ), "data/gtf")
gtf_filename = 'Mus_musculus_casteij.CAST_EiJ_v1.108.gtf'
gff_dir = file.path(getwd( ), "data/gff")
gff_filename = 'Mus_musculus_casteij.CAST_EiJ_v1.108.chr.gff3.gz'
out_dir = file.path(getwd( ), "data/gtf/strip_cast")
out_filename = "strip_cast-Mus_musculus_casteij.tsv"


# Start Log
log <- log_open(paste("strip_cast ", Sys.time(), '.log', sep='')) # Open log
log_print(paste('gtf file: ', file.path(gtf_dir, gtf_filename)))
log_print(paste('gff file: ', file.path(gff_dir, gff_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

# Import gtf data
log_print("Importing gtf data...")
# gtf_cast <- rtracklayer::import('data/gtf/Mus_musculus.GRCm38.87.gtf')
gtf_cast <- rtracklayer::import('data/gtf/Mus_musculus_casteij.CAST_EiJ_v1.108.gtf')
gtf_cast_df <- as.data.frame(gtf_cast, stringsAsFactors = FALSE)
test_cast <- gtf_cast_df[is.na(gtf_cast_df$gene_name) == FALSE & gtf_cast_df$seqnames == 'X',]
cast_filtered <- test_cast[,c(1,2,3,4,7,10,23)]


# First, import the GTF-file that you have also used as input for htseq-count
start_time = Sys.time()
log_print(paste("Making TxDb from gff...", start_time))
txdb_cast <- makeTxDbFromGFF(
    'data/gff/Mus_musculus_casteij.CAST_EiJ_v1.108.chr.gff3.gz',
    # format="gtf"
)
end_time = Sys.time()
log_print(paste('finished.', end_time))
log_print(paste('time of operation... ', end_time - start_time))


# ----------------------------------------------------------------------
# Processing

# then collect the exons per gene id
exons.list.per.gene.cast <- exonsBy(txdb_cast,by="gene")

# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes.cast <- as.data.frame(sum(width(reduce(exons.list.per.gene.cast))))

# try to match exonic lengths to genes on the X using the cast_filtered df and exonic.gene.sizes df.

strip_cast <- cast_filtered[,c(6,7)]    # Pulling out just gene_ids and gene_names from the gtf file
strip_cast <- unique(strip_cast)        # Stripping down this new df to only unique instances

index_in_exonic.gene.sizes.cast <- match(strip_cast[,1], rownames(exonic.gene.sizes.cast))  # find the indexes of the entries in exonic.gene.sizes that matches the genes in the strip_cast dataframe.

strip_cast[,3] <- exonic.gene.sizes.cast[index_in_exonic.gene.sizes.cast,] # Append those indexes onto the strip_cast df and rename cols. 
colnames(strip_cast) <- c('gene_id', 'gene_name', 'exonic_length')


# ----------------------------------------------------------------------
# Save Output

log_print(paste('writing file...', Sys.time()))
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    strip_cast,
    file = file.path(out_dir, out_filename),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t'
)

log_print(paste('End', Sys.time()))
log_close()
