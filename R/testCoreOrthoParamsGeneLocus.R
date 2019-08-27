#' @title Helper function for plot_diverse_homology_thresholds_core_orthologs
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