## Common use functions

#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}

#' Read all the csv files from a directory and append them into a single dataframe
#' 
#' @export
append_many_csv <- function(dir_path, sep='\t', row_names=1) {
    filenames <- list.files(dir_path, full.names=TRUE)
    csv <- lapply(filenames, read.csv, sep=sep, row.names=row_names)
    data <- do.call(rbind, csv)
    return(data)
}

#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
#' index_cols=c('gene_id', 'gene_name', 'chromosome')
#' index_cols=c('count')
#' 
#' @export
join_many_csv <- function(dir_path, index_cols, value_cols, sep=',') {
    filepaths <- list.files(dir_path, full.names=TRUE)
    if (length(filepaths)==0) {
        stop("no files found")
    }
    filenames = c(tools::file_path_sans_ext(basename(filepaths)))
    
    # read dfs and left join on index_cols
    df_list <- lapply(filepaths, read.csv, sep=sep)
    data <- Reduce(
        function(...) merge(..., by=index_cols),
        lapply(df_list, "[", c(index_cols, value_cols))
    )
    
    # rename columns
    colnames(data) = c(
        index_cols,  # index_cols
        as.list(outer(filenames, value_cols, paste, sep='_'))  # prefix value_cols with filename
    )
    return(data)
}