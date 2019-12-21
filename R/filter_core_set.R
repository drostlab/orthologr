#' @title Helper function to extract a core set of orthologous gene loci
#' @description Helper function to extract a core set of orthologous gene loci
#' @param x input data in \code{data.frame} or \code{tibble} format.
#' @param order_species a character vector containing the scientific names of the organisms of interest
#' ordered according to their phylogenetic distance to their reference species.  
filter_core_set <- function(x, order_species) {
        
        subset_species <- sort(as.character(names(table(
                x$subject_species
        ))))
        general_species <- sort(order_species)
        
        if (identical(subset_species, general_species)) {
                return(x)
        } else {
                return(
                        tibble::tibble(
                                query_species = NA,
                                subject_species = NA,
                                query_id = NA,
                                query_gene_locus_id = NA,
                                subject_id = NA,
                                subject_gene_locus_id = NA,
                                dN = NA,
                                dS = NA,
                                dNdS = NA,
                                evalue = NA,
                                bit_score = NA,
                                perc_identity = NA,
                                num_ident_matches = NA,
                                alig_length = NA,
                                mismatches = NA,
                                gap_openings = NA,
                                n_gaps = NA,
                                pos_match = NA,
                                ppos = NA,
                                q_start = NA,
                                q_end = NA,
                                q_len = NA,
                                qcov = NA,
                                qcovhsp = NA,
                                s_start = NA,
                                s_end = NA
                        )
                )

        }
}