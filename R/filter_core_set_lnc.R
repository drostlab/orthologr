#' @title Helper function to extract a core set of orthologous lncRNAs (deprecated)
#' @description This function is deprecated and will be removed in a future release.
#' Please use \code{\link{filter_core_set}} instead, which now handles both
#' protein-coding gene tables and lncRNA maps automatically.
#' @param x input data in \code{data.frame} or \code{tibble} format.
#' @param order_species a character vector containing the scientific names of the organisms of interest
#' ordered according to their phylogenetic distance to their reference species.
filter_core_set_lnc <- function(x, order_species) {
        .Deprecated("filter_core_set")
        filter_core_set(x, order_species)
}
