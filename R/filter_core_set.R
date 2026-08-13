#' @title Helper function to extract a core set of orthologs
#' @description Helper function to extract a core set of orthologs shared across all
#' subject species. Handles both protein-coding gene tables (species column
#' \code{subject_species}) and lncRNA maps (species column \code{species}).
#' The appropriate column is detected automatically from the input.
#' @param x input data in \code{data.frame} or \code{tibble} format.
#' @param order_species a character vector containing the scientific names of the organisms of interest
#' ordered according to their phylogenetic distance to their reference species.
#' @examples \dontrun{
#' # Protein-coding table: all target species present -> returns x unchanged
#' x_complete <- tibble::tibble(
#'   subject_species     = c("Arabidopsis_lyrata", "Brassica_rapa"),
#'   query_id            = c("AT1G01010.1", "AT1G01010.1"),
#'   query_gene_locus_id = c("AT1G01010", "AT1G01010"),
#'   dN = c(0.05, 0.08), dS = c(0.12, 0.15), dNdS = c(0.42, 0.53)
#' )
#' filter_core_set(x_complete, order_species = c("Arabidopsis_lyrata", "Brassica_rapa"))
#'
#' # Protein-coding table: one species missing -> returns a schema-preserving NA row
#' x_incomplete <- tibble::tibble(
#'   subject_species     = "Arabidopsis_lyrata",
#'   query_id            = "AT1G01010.1",
#'   query_gene_locus_id = "AT1G01010",
#'   dN = 0.05, dS = 0.12, dNdS = 0.42
#' )
#' filter_core_set(x_incomplete, order_species = c("Arabidopsis_lyrata", "Brassica_rapa"))
#' }
filter_core_set <- function(x, order_species) {

        species_col <- if ("subject_species" %in% names(x)) "subject_species" else "species"

        subset_species <- sort(as.character(names(table(x[[species_col]]))))
        general_species <- sort(order_species)

        if (identical(subset_species, general_species)) {
                return(x)
        } else {
                # Return a schema-preserving single NA row rather than a hardcoded tibble,
                # so this helper works regardless of whether x is a protein-coding or lncRNA table.
                return(x[NA_integer_, ])
        }
}
