#' @title Retrieve a core set of orthologous lncRNAs (deprecated)
#' @description This function is deprecated and will be removed in a future release.
#' Please use \code{\link{retrieve_core_orthologs}} instead, which now handles both
#' protein-coding gene tables and lncRNA maps automatically.
#' @param lnc_map a \code{lnc_map} that was generated with \code{\link{map_generator_lnc}}.
#' @param species_order a character string specifying species names listed in the order of
#' phylogenetic/taxonomic distance from the query species.
#' @author Hajk-Georg Drost
#' @export

lnc_map_core_orthologs <- function(lnc_map, species_order) {
        .Deprecated("retrieve_core_orthologs")
        retrieve_core_orthologs(ortho_tables = lnc_map, species_order = species_order)
}
