## Functions
## dict
## filter_list_for_match
## items_in_a_not_b


#' Dictionary
#' 
#' @description instantiate a named list. I am aware this is not a real dictionary, this method will be updated in the future.
#'
dict <- function(keys, values) {
    return(setNames(values, keys))
}


#' Find matches based on substring
#' 
#' @description Return elements of a list matching a particular substring
#' 
#' @param items list or vector
#' @param patterns a string or collection of strings
#' @return Returns a list or vector with any matching items
#' 
#' @examples
#' filter_list_for_match(c("a_suffix", "b_suffix", "c"), "suffix")
#' 
#' @export
filter_list_for_match <- function(items, pattern) {
    # filter
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


#' Return items unique to vector a
#' 
#' @description
#' Return all the items found in a and not b. Useful for filtering dataframe columns.
#' 
#' @param a list or vector
#' @param b list or vector
#' @return Returns the filtered list or vector (same type as a)
#' 
#' @references
#' \href{https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list}{Stack Overflow}
#' 
#' @examples
#' items_in_a_not_b(c('a', 'b', 'c', '1', '2'), c('1', '2'))
#' 
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}
