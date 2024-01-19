## Functions
## list_files
## append_many_csv
## join_many_csv
## read_csv_or_tsv


#' List all files with a specific extension
#' 
#' @description
#' This is a thin wrapper around [list.files()].
#' 
#' @references
#' \href{https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching}{StackOverflow post}
#' 
#' @export
list_files <- function(dir_path, recursive = TRUE, full_name=TRUE, ext=NULL) {
    all_files = list.files(dir_path, recursive = recursive, full.name=full_name)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Aggregate csv files by appending them rowwise
#' 
#' @description
#' Read all the csv files from a directory and append them into a single dataframe
#' 
#' @export
append_many_csv <- function(dir_path, sep=',', row_names=NULL) {
    filenames <- list.files(dir_path, full.names=TRUE)
    csv <- lapply(filenames, read.csv, sep=sep, row.names=row_names)
    data <- do.call(rbind, csv)
    return(data)
}


#' Aggregate csv files by joining them column-wise
#' 
#' @description
#' Read all the csv files from a directory and left join them into a single dataframe
#' 
#' @references
#' \href{https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-all_reads-frames}{StackOverflow post}
#' 
#' @export
join_many_csv <- function(paths, index_cols, value_cols, recursive=TRUE) {

    # distinguish if paths is a directory or a list of files
    if (length(paths)==1) {
        if (dir.exists(paths)) {
             paths <- list_files(paths, recursive=recursive)
        }
    } else if (length(paths[file.exists(paths)]) == 0) {
        stop(paste("no files found!"))
    }

    # split into csv and tsv files
    csv_paths = filter_list_for_match(paths, 'csv')
    csv_list <- lapply(csv_paths, read.csv, sep=',')
    tsv_paths = filter_list_for_match(paths, 'tsv')
    tsv_list <- lapply(tsv_paths, read.csv, sep='\t') 
    df_list = c(csv_list, tsv_list)

    filenames = c(tools::file_path_sans_ext(basename(paths)))
    # Warning: column names ‘count.x’, ‘count.y’ are duplicated in the result
    # See: https://stackoverflow.com/questions/38603668/suppress-any-emission-of-a-particular-warning-message
    withCallingHandlers({
        merged <- Reduce(
            function(...) merge(..., by=index_cols),
            lapply(df_list, "[", c(index_cols, value_cols))
        )
    }, warning = function(w) {
        # print(conditionMessage(w))
        if (startsWith(conditionMessage(w), "column names")) {
            invokeRestart( "muffleWarning" )
        }
    })
    
    # rename columns
    colnames(merged) = c(
        index_cols,  # index_cols
        as.list(outer(value_cols, filenames, paste, sep='-'))  # suffix value_cols with filename
    )
    return(merged)
}


#' Switch case to read csv or tsv based on the extension
#'
#' @description Mainly used to simplify scripts
#' 
#' @export
read_csv_or_tsv <- function(file, header=TRUE, sep=',', check_names=FALSE) {
    ext = tools::file_ext(file)
    if (ext == 'tsv') {
        sep='\t'
    }

    data = read.csv(file, header=header, sep=sep, check.names=check_names)
    return(data)
}
