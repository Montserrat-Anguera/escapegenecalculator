import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'join_many_csv', 'read_csv_or_tsv', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)

## Functions
## get_single_reads
## get_merged_reads
## concat_reads


#' Convenience function
#' 
get_single_reads <- function(
    reads_filename,
    chromosomal_parentage='mat',
    mouse_strain='Mus_musculus'
) {
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


#' Convenience function to only read relevant rows into memory
#'
get_merged_reads <- function(reads_dir, config, chromosome_parentage='mat') {

    index_cols_ = c('gene_id', 'gene_name', 'chromosome')
    value_cols_ = 'count'
    reads <- join_many_csv(reads_dir, index_cols=index_cols_, value_cols=value_cols_)

    colnames(reads) = c(
        index_cols_,
        paste0('count', '-', chromosome_parentage, '-', config[['mouse_id']])
    )

    # standardize columns
    colnames(reads) <- tolower(colnames(reads))
    colnames(reads) <- lapply(
      colnames(reads),
      function(x) {
        step1 <- gsub('-', '_', x)  # convert mixed_case-col_names to fully snake_case
        step2 <- gsub('^chr_', 'chromosome_', step1)
        step3 <- gsub('^count_', 'num_reads_', step2)
        step4 <- gsub('^reads_', 'num_reads_', step3)
        return (step4)}
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

    # example dataframe
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
