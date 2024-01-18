## Functions
## camel_to_snake_case


#' Standardize Camel-case Titles
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