#' @title Helper function to select best BLAST hit based on min evalue
filter_best_hits <- function(x) {
        evalue <- NULL
        res <- dplyr::filter(x, evalue == min(evalue))
        return(res)
}