import::here(tidyr, 'pivot_wider')

## Functions
## filter_dataframe_column_by_list
## most_frequent_item
## pivot
## reset_index


#' Filter dataframe column by list
#' 
#' @description
#' Use this to search rows for a list of values.
#' 
#' @param df a dataframe
#' @param colname select a column
#' @param items used to filter
#' @return Returns a dataframe with rows in which items match the column
#' 
#' @examples
#' filter_dataframe_column_by_list(mtcars, 'cyl',  c(4, 6))
#' 
#' @export
filter_dataframe_column_by_list <- function(dataframe, colname, items, index_name='index', return_index=FALSE) {
    
    data <- reset_index(dataframe, index_name=index_name)
    rownames(data) <- data[, colname]
    data <- (data[intersect(data[, colname], items),])  # filter genes shared by both gtf files
    rownames(data) <- data[, index_name]  # optional preserve index for troubleshooting
    
    if (return_index==TRUE) {
        return (data)
    } else {
        return (data[, items_in_a_not_b(colnames(data), 'index')])
    }
}


#' Get the most frequently occurring item in a dataframe column
#'
#' @export
most_frequent_item <- function(df, colname) {
    return(df[which.max(factor(df[, colname])), colname])
}


#' Pivot
#' 
#' @description Thin wrapper around [tidyr::pivot_wider()]. Return a native dataframe instead of a tibble
#' 
#' @examples
#' pivot(reset_index(mtcars), columns='index', values='mpg')
#' 
#' @seealso
#' [tidyr::pivot_wider()]
#' 
#' @export
pivot <- function(df, columns, values) {

    # Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    # See: http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings
    withCallingHandlers({
        tibble_obj = tidyr::pivot_wider(
            df,
            names_from = columns,
            values_from = values,
            values_fn = list,
            names_glue = if (length(values)==1){"{.value}_{.name}"}
        )
    }, warning = function(w) {
        if ( any(grepl("Using an external vector", w)) ) {
            invokeRestart("muffleWarning")
        }
    })

    # Note: This duplicates the functionality of data.table::setDT(tibble_obj)
    # Note 2: The reverse operation is tibble::as_tibble(tibble_obj)
    dataframe = data.frame(lapply(tibble_obj, unlist))
    
    return (dataframe)
}


#' Reset index
#' 
#' @description
#' Moves the values in the index to a column. Resets the index to the default integer index.
#' Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.reset_index.html}{reset_index}.
#' 
#' @param df a dataframe
#' @param index_name select a new name for the index column
#' @param drop if TRUE, does not copy the index values to the new column
#' @return Returns a dataframe with index renamed and new integer index 
#' 
#' @examples
#' reset_index(mtcars, index_name='model')
#' 
#' @references
#' \href{https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column}{StackOverflow post}
#' 
#' @exportå
reset_index <- function(df, index_name='index') {
    df <- cbind(index = rownames(df), df)
    rownames(df) <- 1:nrow(df)
    colnames(df)[colnames(df) == "index"] = index_name
    return (df)
}
