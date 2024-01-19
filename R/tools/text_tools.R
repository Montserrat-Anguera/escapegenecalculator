import::here(stringr, 'str_to_title')

## Functions
## camel_to_snake_case
## snake_to_title_case


#' Standardize CamelCase Titles
#' 
#' @description Converts "ColumnTitle" to column_title
#' 
#' @examples
#' camel_to_snake_case('ColumnTitle')
#' 
#' @export
camel_to_snake_case <- function(text) {
    return(tolower(gsub("([a-z])([A-Z])", "\\1_\\L\\2", text, perl = TRUE)))
}


#' Standardize snake_case Titles
#' 
#' @description Converts column_title to "Column Title"
#' 
#' @examples
#' snake_to_title_case('column_title')
#' 
#' @export
snake_to_title_case <- function(text) {
    return(str_to_title(gsub("_", " ", text)))
}
