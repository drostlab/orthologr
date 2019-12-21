#' @title Helper function for \code{plot_diverse_homology_thresholds_core_orthologs}
#' @description Helper function for \code{\link{plot_diverse_homology_thresholds_core_orthologs}}.
#' @param ortho_tables a \code{ortho tables} that was generated with \code{\link{generate_ortholog_tables_all}}.
#' @param param parameters.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
testCoreOrthoParamsGeneLocus <-
        function(ortho_tables, param, species_order) {
                ortho_tables_core <-
                        orthologr::retrieve_core_orthologs(ortho_tables = ortho_tables,
                                                           species_order = species_order)
                
                query_gene_locus_id <- NULL
                return(tibble::tibble(
                        parameter = param,
                        n_core_genes = nrow(
                                dplyr::summarize(
                                        dplyr::group_by(ortho_tables_core, query_gene_locus_id),
                                        n_genes = dplyr::n()
                                )
                        )
                ))
        }