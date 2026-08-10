#' @title Retrieve a core set of orthologs from pairwise ortholog tables
#' @description Given an ortholog table generated with \code{\link{generate_ortholog_tables_all}}
#' or a lncRNA map generated with \code{\link{map_generator_lnc}}, this function determines
#' a core set of orthologs that are shared between all species.
#'
#' Both table types are supported and detected automatically:
#' \itemize{
#'   \item Protein-coding tables use species column \code{subject_species} and query column
#'         \code{query_gene_locus_id}.
#'   \item lncRNA maps use species column \code{species} and query column \code{query_id}.
#' }
#' @param ortho_tables an ortholog table generated with \code{\link{generate_ortholog_tables_all}}
#' or a lncRNA map generated with \code{\link{map_generator_lnc}}.
#' @param species_order a character string specifying species names listed in the order of
#' phylogenetic/taxonomic distance from the query species. The species names must match
#' the species names present in \code{ortho_tables}.
#' @author Hajk-Georg Drost
#' @examples
#' # Protein-coding ortholog table (subject_species + query_gene_locus_id)
#' ortho_tbl <- tibble::tibble(
#'   query_species       = "Arabidopsis_thaliana",
#'   subject_species     = rep(c("Arabidopsis_lyrata", "Brassica_rapa"), each = 3),
#'   query_id            = rep(c("AT1G01010.1", "AT1G01020.1", "AT1G01030.1"), 2),
#'   query_gene_locus_id = rep(c("AT1G01010", "AT1G01020", "AT1G01030"), 2),
#'   subject_id          = paste0("subj_", seq_len(6)),
#'   q_len               = rep(c(430L, 246L, 359L), 2),
#'   alig_length         = rep(c(430L, 246L, 355L), 2)
#' )
#' retrieve_core_orthologs(ortho_tbl,
#'                         species_order = c("Arabidopsis_lyrata", "Brassica_rapa"))
#'
#' # lncRNA map (species + query_id) — column layout detected automatically
#' lnc_tbl <- tibble::tibble(
#'   species     = rep(c("Arabidopsis_lyrata", "Brassica_rapa"), each = 2),
#'   query_id    = rep(c("lnc001", "lnc002"), 2),
#'   subject_id  = paste0("slnc_", seq_len(4)),
#'   q_len       = rep(c(500L, 800L), 2),
#'   alig_length = rep(c(490L, 800L), 2)
#' )
#' retrieve_core_orthologs(lnc_tbl,
#'                         species_order = c("Arabidopsis_lyrata", "Brassica_rapa"))
#' @export

retrieve_core_orthologs <- function(ortho_tables, species_order) {

        species_col <- if ("subject_species" %in% names(ortho_tables)) "subject_species" else "species"
        query_col   <- if ("query_gene_locus_id" %in% names(ortho_tables)) "query_gene_locus_id" else "query_id"

        if (length(names(table(ortho_tables[[species_col]]))) != length(species_order))
                stop("The number of different species specified in 'ortho_tables' does not match with the number of different species specified in 'species_order'.", call. = FALSE)

        message("Retrieving core orthologs that are present in all pairwise species comparisons: ",
                paste0(names(table(ortho_tables[[species_col]])), collapse = ", "))

        q_len <- alig_length <- NULL
        ortho_tables <- dplyr::mutate(
                ortho_tables,
                scope = 1 - (abs(q_len - alig_length) / q_len)
        )

        query_species_label <- if ("query_species" %in% names(ortho_tables))
                paste0("from ", unique(ortho_tables$query_species), " ")
        else
                ""

        message("The species order in terms of phylogenetic distance ", query_species_label,
                "was set to: ", paste0(species_order, collapse = ", "))

        ortho_table_core <- dplyr::bind_rows(dplyr::group_map(
                .data = dplyr::group_by(ortho_tables, .data[[query_col]]),
                .f = ~ filter_core_set(., order_species = species_order),
                .keep = TRUE
        ))

        ortho_table_core <- dplyr::filter(ortho_table_core, !is.na(.data[[query_col]]))

        core_grouped <- dplyr::summarize(
                dplyr::group_by(ortho_table_core, .data[[query_col]]),
                n = dplyr::n()
        )

        if (!all(core_grouped$n == length(species_order)))
                stop(paste0("Somehow there seem to be non-unique ", query_col,
                            "'s in the core set. Please check what could have gone wrong."),
                     call. = FALSE)

        message("Core orthologs (n = ", nrow(core_grouped), " loci) were successfully retrieved.")
        return(ortho_table_core)
}
