read.ortholog.table <- function(file) {
        
        res <- readr::read_delim(
                file,
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
                        "alig_length" = readr::col_integer(),
                        "mismatches" = readr::col_integer(),
                        "gap_openings" = readr::col_integer(),
                        "q_start" = readr::col_integer(),
                        "q_end" = readr::col_integer(),
                        "s_start" = readr::col_integer(),
                        "s_end" = readr::col_integer()
                )
        )
        
        return(res)
}