#' @title Importing output pairwise orthologs tables generated with \code{\link{generate_ortholog_tables_all}}
#' @description Importing output pairwise orthologs tables generated with \code{\link{generate_ortholog_tables_all}}.
#' @param ortholog_tables_folder file path to output files generated with \code{\link{generate_ortholog_tables_all}}.
#' @author Hajk-Georg Drost
#' @export
import_ortholog_tables_all <- function(ortholog_tables_folder) {
        
        if (!file.exists(ortholog_tables_folder))
                stop("The folder '",ortholog_tables_folder, "' does not seem to exist. Please provide a valid folder path.", call. = FALSE)
        
        message("Importing ortholog tab les generated with generate_ortholog_tables_all() from '", ortholog_tables_folder, " ...")
        
        ortholog_tables <- file.path(ortholog_tables_folder, list.files(ortholog_tables_folder))
        
        res <-
                suppressMessages(dplyr::bind_rows(lapply(ortholog_tables, function(x)
                        readr::read_delim(
                                x,
                                col_names = TRUE,
                                delim = ";",
                                col_types = readr::cols(
                                        "query_species" = readr::col_character(),
                                        "subject_species" = readr::col_character(),
                                        "query_id" = readr::col_character(),
                                        "query_gene_locus_id" = readr::col_character(),
                                        "subject_id" = readr::col_character(),
                                        "subject_gene_locus_id" = readr::col_character(),
                                        "dN" = readr::col_double(),
                                        "dS"  = readr::col_double(),
                                        "dNdS" = readr::col_double(),
                                        "evalue" = readr::col_double(),
                                        "bit_score" = readr::col_double(),
                                        "perc_identity" = readr::col_double(),
                                        "num_ident_matches" = readr::col_double(),
                                        "alig_length" = readr::col_integer(),
                                        "mismatches" = readr::col_integer(),
                                        "gap_openings" = readr::col_integer(),
                                        "n_gaps" = readr::col_double(),
                                        "pos_match" = readr::col_double(),
                                        "ppos" = readr::col_double(),
                                        "q_start" = readr::col_integer(),
                                        "q_end" = readr::col_integer(),
                                        "q_len" = readr::col_double(),
                                        "qcov" = readr::col_double(),
                                        "qcovhsp" = readr::col_double(),
                                        "s_start" = readr::col_integer(),
                                        "s_end" = readr::col_integer(),
                                        "s_len" = readr::col_double(),
                                        "score_raw" = readr::col_double()
                                )
                        ))))
        message("Import was successful!")
        return(res)
}