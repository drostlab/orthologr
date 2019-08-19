#' @title Retrieve a core set of orthologous lncRNAs from the pairwise lncRNA orthologs map
#' @description Given a lnc_map table generated with \code{\link{map_generator_lnc}},
#' this function will determine a core set of lncRNA orthologs that are shared between all species.
#' @param ortho_tables a \code{lnc_map} that was generated with \code{\link{map_generator_lnc}}.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in the \code{lnc_map}.
#' @author Hajk-Georg Drost
#' @export

lnc_map_core_orthologs <- function(lnc_map, species_order) {
        
        if (length(names(table(lnc_map$species))) != length(species_order))
                stop("The number of different species specified in 'lnc_map' does not match with the number of different species specified in 'species_order'.", call. = FALSE)
        
        message("Retrieving core lncRNA orthologs that are present in all pairwise species comparisons: ", paste0(names(table(lnc_map$species)), collapse = ", "))
        q_len <- alig_length <- query_id <- NULL
        lnc_map <-
                dplyr::mutate(
                        lnc_map,
                        scope = 1 - (abs(q_len - alig_length) / q_len)
                )
        message("The species order in terms of phylogenetic distance from the query species was set to: ", paste0(species_order, collapse = ", "))
        
        lnc_map_core <-
                dplyr::bind_rows(dplyr::group_map(
                        .tbl = dplyr::group_by(lnc_map, query_id),
                        ~ filter_core_set_lnc(., order_species = species_order),
                        keep = TRUE
                ))
        
        lnc_map_core <-
                dplyr::filter(lnc_map_core, !is.na(query_id))
        
        core_grouped <- dplyr::summarize(dplyr::group_by(lnc_map_core, query_id), n = dplyr::n())
        
        if (!all(core_grouped$n == length(species_order)))
                stop("Somehow there seem to be non-unique query_id's in the core set. Please check what could have gone wrong.", call. = FALSE)
        
        message("Core orthologs (n = ", nrow(core_grouped), " loci) were successfully retrieved.")
        return(lnc_map_core)
}
