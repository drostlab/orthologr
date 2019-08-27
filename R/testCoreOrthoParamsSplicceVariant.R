#' @title Helper function for plot_diverse_homology_thresholds_core_orthologs 2
testCoreOrthoParamsSpliceVariant <-
        function(ortho_tables,
                 param,
                 species_order,
                 n_core_species) {
                ortho_tables_core <-
                        orthologr::retrieve_core_orthologs(ortho_tables = ortho_tables,
                                                           species_order = species_order)
                
                query_id <- n_genes <- NULL
                return(tibble::tibble(
                        parameter = param,
                        n_core_genes = nrow(
                                dplyr::filter(
                                        dplyr::summarize(
                                                dplyr::group_by(ortho_tables_core, query_id),
                                                n_genes = dplyr::n()
                                        ),
                                        n_genes == n_core_species
                                )
                        )
                ))
        }
