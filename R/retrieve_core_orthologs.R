#' @title Retrieve a core set of orthologous gene loci from several pairwise ortholog tables
#' @description Given a ortho table generated with \code{\link{generate_ortholog_tables_all}},
#' this function will determine a core set of orthologs that are shared between all species.
#' @param ortho_tables a \code{ortho tables} that was generated with \code{\link{generate_ortholog_tables_all}}.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in the \code{ortho_tables}.
#' @author Hajk-Georg Drost
#' @export

retrieve_core_orthologs <- function(ortho_tables, species_order) {
        
        if (length(names(table(ortho_tables$subject_species))) != length(species_order))
                stop("The number of different species specified in 'ortho_tables' does not match with the number of different species specified in 'species_order'.", call. = FALSE)
        
        message("Retrieving core orthologs that are present in all pairwise species comparisons: ", paste0(names(table(ortho_tables$subject_species)), collapse = ", "))
        q_len <- alig_length <- query_gene_locus_id <- NULL
        ortho_tables <-
                dplyr::mutate(
                        ortho_tables,
                        scope = 1 - (abs(q_len - alig_length) / q_len)
                )
        message("The species order in terms of phylogenetic distance from ",unique(ortho_tables$query_species)," was set to: ", paste0(species_order, collapse = ", "))
     
        ortho_table_core <-
                dplyr::bind_rows(dplyr::group_map(
                        .data = dplyr::group_by(ortho_tables, query_gene_locus_id),
                        .f = ~ filter_core_set(., order_species = species_order),
                        .keep = TRUE
                ))
        
        ortho_table_core <-
                dplyr::filter(ortho_table_core, !is.na(query_gene_locus_id))

        core_grouped <- dplyr::summarize(dplyr::group_by(ortho_table_core, query_gene_locus_id), n = dplyr::n())
        
        if (!all(core_grouped$n == length(species_order)))
                stop("Somehow there seem to be non-unique query_gene_locus_id's in the core set. Please check what could have gone wrong.", call. = FALSE)

        message("Core orthologs (n = ", nrow(core_grouped), " loci) were successfully retrieved.")
        return(ortho_table_core)
}
