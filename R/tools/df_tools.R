import::here(tidyr, 'pivot_wider')

## Functions
## coalesce1
## coalesce_colnames
## fillna
## filter_dataframe_column_by_list
## multi_melt
## most_frequent_item
## pivot
## reset_index


#' Coalesce
#' 
#' @references
#' \href{https://stackoverflow.com/questions/19253820/how-to-implement-coalesce-efficiently-in-r}{StackOverflow post}
#' 
#' 
coalesce1 <- function(...) {
    ans <- ..1
    for (elt in list(...)[-1]) {
        i <- is.na(ans)
        ans[i] <- elt[i]
    }
    ans
}


#' Coalesce Column Names
#' 
#' @description
#' Given a dataframe with columns with one-hot encodings,
#' aggregates the column names into a comma-separated list
#' 
#' @export
coalesce_colnames <- function(df, cols, sep=', ') {

    output <- ''
    for (col in cols) {
        output <- gsub('1', col, paste(output, gsub('0', '', df[[col]]), sep=sep))
        output <- gsub(' ,', '', output)
    }
    output <- sub('^, ', '', output)
    output <- sub(', $', '', output)
    return(output)
}


#' Fill specific column with NA
#' 
#' @description Mirrors Pandas' \href{https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.fillna.html}{fillna}
#' 
#' @param df a dataframe
#' @param cols a list of columns to replace
#' @param val the value to fill with
#' @param inplace TRUE allows you to avoid re-assigning the variable
#' @return Returns a dataframe.
#' 
#' @examples
#' mtcars['new_col'] <- NA
#' head(fillna(mtcars, c('new_col'), 1))
#' 
#' @export
fillna <- function(df, cols, val=0, inplace=FALSE) {
    df_name <- deparse(substitute(df))
    for (col in cols) {
        df[is.na(df[, col]), col] <- val
    }
    if (inplace) {
        assign(df_name, df, envir=.GlobalEnv)
    } else {
        return(df)
    }
}


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


#' Multiple melt
#' 
#' @description Unpivot multiple columns simulataneously
#' 
#' @export
multi_melt <- function(
    df,
    id_vars=c('id', 'name'),
    values=list(
        'metric_1'=c('metric_1_x', 'metric_1_y'),
        'metric_2'=c('metric_2_x', 'metric_2_y')
    ),
    var_name='variable'
) {

    dfs = new.env()
    for (value_name in names(values)) {
        value_vars <- values[[value_name]]
        tmp = reshape2::melt(
            df[, c(id_vars, value_vars)],
            row=id_vars,
            measure.vars=value_vars,
            value.name=value_name,
            variable.name=var_name
        )

        tmp[var_name] = sub(paste0(value_name, '_'), '', tmp[[var_name]])
        
        dfs[[value_name]] <- tmp
    }
    dfs <- as.list(dfs)

    merged <- Reduce(
        function(...) merge(..., by=(c(id_vars, var_name))),
        dfs
    )

    return(merged)
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
