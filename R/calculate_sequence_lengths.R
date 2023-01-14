# 
# # Import the gzipped gtf file using rtracklayer, convert it to a data frame, throw away rows that have an NA in the gene_name column, throw out rows that have genes mapped to autosomes, and pull out the chromosome, start, stop, width, type, gene id, and gene name columns. Then loop through and pull out the width and gene names that correspond to type == 'gene'. Do for both genomes. 
# 
# gtf_cast <- rtracklayer::import('/Users/sarahpyfrom/Box\ Sync/ANGUERA/Genome_Files/GTF/Mus_musculus.GRCm38.87.gtf')
# gtf_cast_df <- as.data.frame(gtf_cast, stringsAsFactors = FALSE)
# test_cast <- gtf_cast_df[is.na(gtf_cast_df$gene_name) == FALSE & gtf_cast_df$seqnames == 'X',]
# cast_filtered <- test_cast[,c(1,2,3,4,7,10,23)]
# 
# 
# # First, import the GTF-file that you have also used as input for htseq-count
# 
# library(GenomicFeatures)
# txdb_cast <- makeTxDbFromGFF('Mus_musculus_casteij.CAST_EiJ_v1.97.gtf.gz',format="gtf")
# 
# # then collect the exons per gene id
# 
# exons.list.per.gene.cast <- exonsBy(txdb_cast,by="gene")
# 
# # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
# 
# exonic.gene.sizes.cast <- as.data.frame(sum(width(reduce(exons.list.per.gene.cast))))
# 
# # try to match exonic lengths to genes on the X using the cast_filtered df and exonic.gene.sizes df.
# 
# strip_cast <- cast_filtered[,c(6,7)]    # Pulling out just gene_ids and gene_names from the gtf file
# strip_cast <- unique(strip_cast)        # Stripping down this new df to only unique instances
# 
# index_in_exonic.gene.sizes.cast <- match(strip_cast[,1], rownames(exonic.gene.sizes.cast))  # find the indexes of the entries in exonic.gene.sizes that matches the genes in the strip_cast dataframe.
# 
# strip_cast[,3] <- exonic.gene.sizes.cast[index_in_exonic.gene.sizes.cast,] # Append those indexes onto the strip_cast df and rename cols. 
# colnames(strip_cast) <- c('gene_id', 'gene_name', 'exonic_length')
# 

# Do the same thing as above for the 129 genome. 
# 
# gtf_b6 <- rtracklayer::import('Mus_musculus_c57bl6nj.C57BL_6NJ_v1.96.gtf.gz')
# gtf_b6_df <- as.data.frame(gtf_b6, stringsAsFactors = FALSE)
# test_b6 <- gtf_b6_df[is.na(gtf_b6_df$gene_name) == FALSE & gtf_b6_df$seqnames == 'X',]
# b6_filtered <- test_b6[,c(1,2,3,4,7,10,23)]
# 
# # First, import the GTF-file that you have also used as input for htseq-count
# 
# #library(GenomicFeatures)
# txdb_b6 <- makeTxDbFromGFF('Mus_musculus_c57bl6nj.C57BL_6NJ_v1.96.gtf.gz',format="gtf")
# 
# # then collect the exons per gene id
# 
# exons.list.per.gene.b6 <- exonsBy(txdb_b6,by="gene")
# 
# # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
# 
# exonic.gene.sizes.b6 <- as.data.frame(sum(width(reduce(exons.list.per.gene.b6))))
# 
# # try to match exonic lengths to genes on the X using the cast_filtered df and exonic.gene.sizes df.
# 
# strip_b6 <- b6_filtered[,c(6,7)]    # Pulling out just gene_ids and gene_names from the gtf file
# strip_b6 <- unique(strip_b6)        # Stripping down this new df to only unique instances
# 
# index_in_exonic.gene.sizes.b6 <- match(strip_b6[,1], rownames(exonic.gene.sizes.b6))  # find the indexes of the entries in exonic.gene.sizes that matches the genes in the strip_cast dataframe.
# 
# strip_b6[,3] <- exonic.gene.sizes.b6[index_in_exonic.gene.sizes.b6,] # Append those indexes onto the strip_cast df and rename cols. 
# colnames(strip_b6) <- c('gene_id', 'gene_name', 'exonic_length')



#Rewriting the above to more directly calculate sequence lengths

library(seqinr)

print('reading fasta file')
fasta_file <- read.fasta(
    file.path(getwd( ), "data/transcriptome", "Mus_musculus_casteij.CAST_EiJ_v1.dna_rm.toplevel.fa")
)

fasta_length <- getLength(fasta_file)
names(fasta_length) <- names(fasta_file)
fasta_length <- data.frame(fasta_length)

print('writing fasta lengths')
write.table(
	rpkms_and_srpms_filtered,
	file = file.path(getwd( ), "data/transcriptome/length", "Mus_musculus_casteij_lengths.tsv")
	row.names = FALSE,
	sep = '\t'
)
