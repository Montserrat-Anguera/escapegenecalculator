
# Import the gzipped gtf file using rtracklayer, convert it to a data frame, throw away rows that have an NA in the gene_name column, throw out rows that have genes mapped to autosomes, and pull out the chromosome, start, stop, width, type, gene id, and gene name columns. Then loop through and pull out the width and gene names that correspond to type == 'gene'. Do for both genomes. 

gtf_cast <- rtracklayer::import('data/gtf/Mus_musculus.GRCm38.87.gtf')
gtf_cast_df <- as.data.frame(gtf_cast, stringsAsFactors = FALSE)
test_cast <- gtf_cast_df[is.na(gtf_cast_df$gene_name) == FALSE & gtf_cast_df$seqnames == 'X',]
cast_filtered <- test_cast[,c(1,2,3,4,7,10,23)]


# First, import the GTF-file that you have also used as input for htseq-count

library(GenomicFeatures)
txdb_cast <- makeTxDbFromGFF('Mus_musculus_casteij.CAST_EiJ_v1.97.gtf.gz',format="gtf")

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
