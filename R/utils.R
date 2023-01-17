## Common use functions

#' Read all the csv files from a directory and append them into a single dataframe
#' 
#' @export
read_many_csv <- function(dir_path, sep='\t', row_names=1) {
    filenames <- list.files(dir_path, full.names=TRUE)
    csv <- lapply(filenames, read.csv, sep=sep, row.names=row_names)
    data <- do.call(rbind, csv)
    return(data)
}

#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
#' by='id' or by=c('gene_id', 'gene_name', 'chromosome')
#' 
#' @export
join_many_csv <- function(dir_path, by, sep='\t') {
    filenames <- list.files(dir_path, full.names=TRUE)
    df_list <- lapply(filenames, read.csv, sep=sep)
    data <- Reduce(function(...) merge(..., by=by), df_list)
    colnames(data) = c(by, c(tools::file_path_sans_ext(basename(filenames))))  # should append
    return(data)
}


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}
