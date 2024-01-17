## Functions
## get_reads
## concat_reads


#' Convenience function
#' 
get_reads <- function(
    reads_filename,
    chromosomal_parentage='mat',
    mouse_strain='Mus_musculus'
) {
    if (!is.null(mouse_strain)) { mouse_strain="Mus_musculus" }  # default option

    reads_file = file.path(in_dir, paste(chromosomal_parentage, "_reads", sep=''), reads_filename)
    reads = read_csv_or_tsv(reads_file)

    exon_lengths_file = file.path(wd, "ref", paste("exon_lengths-", mouse_strain, ".csv", sep=''))
    if (!file.exists(exon_lengths_file)) { return(reads) }

    exon_lengths = read.csv(exon_lengths_file, header=TRUE, check.names=FALSE)
    reads = merge(
        reads,
        exon_lengths[c('gene_name', 'exon_length')],  # do we need 'gene_id'?
        by='gene_name', suffixes=c('', '_'),
        all.x=TRUE, all.y=FALSE
        # na_matches = 'never'  # do not include null values
    )
    return(reads)
}


#' Convenience function to left join summed reads
#'
concat_reads <- function(
        df, reads, cols,
        new_colname="reads", chromosome_parentage="mat", prefix='num_reads_'
    ) {

    # sum the reads
    reads_for_mouse_id = colSums(reads[cols])

    # remove 'count-mat-' from 'count-mat-mouse_1'
    names(reads_for_mouse_id) = gsub(
        paste(prefix, chromosome_parentage, '_', sep=''), '',
        names(reads_for_mouse_id)
    )

    # instantiate the following dataframe
    # +----------+-----------------+
    # | mouse_id | total_reads_mat |
    # +----------+-----------------+
    # | mouse_1  |      4773667    |
    # | mouse_2  |      5222637    |
    # | mouse_3  |      3279278    |
    # +----------+-----------------+
    reads_df = data.frame(
        mouse_id = names(reads_for_mouse_id),
        reads = reads_for_mouse_id
    )
    colnames(reads_df) = c('mouse_id', new_colname)  # rename new column

    # left join to the main table
    df = merge(
        df, reads_df,
        by='mouse_id',
        all.x=FALSE, all.y=FALSE,  # do not include null values
        na_matches = 'never'
    )
    return(df)
}

