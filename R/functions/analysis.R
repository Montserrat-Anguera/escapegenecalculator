## Functions
## adjust_closer


#' Adjust closer
#' 
#' @description
#' If exon_length_pat is too far off from exon_length_mat, adjust it closer
#' This is just for visual aid
#' 
adjust_closer <- function(comparator, base_value, tolerance=0.063) {

    # ignore NA
    if (is.na(comparator) | is.na(base_value)) {
        return(NA)
    }

    if (comparator / base_value <= 1-tolerance) {
        return(base_value * (1-tolerance))
    }
    else if (comparator / base_value >= 1+tolerance) {
        return(base_value * (1+tolerance))
    }
    return(comparator)
}
