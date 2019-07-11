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
        
        res <- suppressMessages(dplyr::bind_rows(lapply(ortholog_tables, function(x) readr::read_delim(x, col_names = TRUE, delim = ";"))))
        message("Import was successful!")
        return(res)
}