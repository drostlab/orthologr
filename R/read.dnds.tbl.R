#' @title Import a dNdS table generated with \code{dNdS}
#' @description This function reads a file that stores a dnds table
#' generated with \code{\link{dNdS}}.
#' @param file file path to dNdS table with \code{;} separated columns.
#' @author Hajk-Georg Drost
#' @examples 
#' # generate dNdS table
#' dNdS_tbl <- dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#' subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#' ortho_detection = "RBH", 
#' aa_aln_type     = "pairwise",
#' aa_aln_tool     = "NW", 
#' codon_aln_tool  = "pal2nal", 
#' dnds_est.method = "Comeron", 
#' comp_cores      = 1 )
#' 
#' # save dNdS table as ';' column separated file
#' utils::write.table(
#' dNdS_tbl, 
#' file.path(tempdir(), "dNdS_tbl.csv"), sep = ";", 
#' col.names = TRUE,
#' row.names = FALSE,
#' quote     = FALSE )
#' 
#' # import dNdS table into R session
#' dNdS_tbl_import <- read.dnds.tbl(file.path(tempdir(), "dNdS_tbl.csv"))
#' @export
read.dnds.tbl <- function(file) {
        
        if (!file.exists(file))
                stop("The file '",file,"' does not seem to exist. Please specify a valid path to your dnds table.", call. = FALSE)
        
        res <- tibble::as_tibble(readr::read_delim(
                file,
                col_names = TRUE,
                delim = ";",
                col_types = readr::cols("query_id" = readr::col_character(),
                                        "subject_id"= readr::col_character(),
                                        "dN" = readr::col_double(),
                                        "dS"  = readr::col_double(),
                                        "dNdS" = readr::col_double(),
                                        "perc_identity" = readr::col_double(),
                                        "alig_length" = readr::col_integer(),
                                        "mismatches" = readr::col_integer(),
                                        "gap_openings" = readr::col_integer(),
                                        "q_start" = readr::col_integer(),
                                        "q_end" = readr::col_integer(),
                                        "s_start" = readr::col_integer(),
                                        "s_end" = readr::col_integer(),
                                        "evalue" = readr::col_double(),
                                        "bit_score" = readr::col_double()
                                        ))
        )
        
        return(res)
}