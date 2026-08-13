#' @title Main Orthology Inference Function
#' @description This function takes nucleotide or protein sequences for a set of organisms
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' @param outgroup_file currently unused; reserved for future use.
#' @param eval a numeric value specifying the E-Value cutoff for hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are: "cds", "protein", or "dna". In case of "cds", sequences are translated to protein sequences,
#' in case of "dna", CDS prediction is performed on the corresponding sequences which are subsequently
#' translated to protein sequences. Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta". Default is "fasta".
#' @param ortho_detection a character string specifying the orthology inference method to use.
#' Default is \code{ortho_detection} = \code{"DIAMOND_RBH"} (DIAMOND2 reciprocal best hit).
#' Available methods:
#' \itemize{
#'   \item \code{"DIAMOND_RBH"}: DIAMOND2 reciprocal best hit (default; fast, recommended)
#'   \item \code{"DIAMOND_BH"}: DIAMOND2 best hit
#'   \item \code{"RBH"}: BLAST reciprocal best hit
#'   \item \code{"BH"}: BLAST best hit
#'   \item \code{"Orthofinder2"}: OrthoFinder2 (currently under development)
#' }
#' @param sensitivity_mode a character string specifying the DIAMOND2 sensitivity level.
#' Only used when \code{ortho_detection} is \code{"DIAMOND_RBH"} or \code{"DIAMOND_BH"}.
#' Options: \code{"faster"}, \code{"fast"} (default), \code{"mid-sensitive"},
#' \code{"sensitive"}, \code{"more-sensitive"}, \code{"very-sensitive"}, \code{"ultra-sensitive"}.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep (DIAMOND2 only).
#' Default is \code{10000}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should
#' be removed from the input \code{file}. Default is \code{FALSE}.
#' @param cdd.path currently unused; reserved for future use.
#' @param path a character string specifying the path to the orthology inference tool executable.
#' @param add_params a character string specifying additional parameters passed to the orthology inference tool.
#' For DIAMOND2 methods, use \code{diamond_params} instead. Default is \code{NULL}.
#' @param diamond_params a character string specifying additional parameters passed to DIAMOND2.
#' Only used when \code{ortho_detection} is \code{"DIAMOND_RBH"} or \code{"DIAMOND_BH"}.
#' Default is \code{NULL}.
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore computations.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param clean_folders a boolean value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @details
#' This function takes sequence files of a query organism and a subject organism and performs orthology inference
#' using a defined orthology inference method to detect orthologous genes.
#'
#' The following interfaces are implemented in the \code{orthologs} function:
#'
#' DIAMOND2 based methods (recommended):
#' \itemize{
#'   \item DIAMOND2 reciprocal best hit (\code{"DIAMOND_RBH"}) — default
#'   \item DIAMOND2 best hit (\code{"DIAMOND_BH"})
#' }
#'
#' BLAST based methods:
#' \itemize{
#'   \item BLAST reciprocal best hit (\code{"RBH"})
#'   \item BLAST best hit (\code{"BH"})
#' }
#'
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of
#' orthologous genes in the second column and the amino acid sequences in the third column.
#' @examples \dontrun{
#'
#' ### DIAMOND2 Reciprocal Best Hit (default)
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein")
#'
#'
#' ### DIAMOND2 Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "DIAMOND_BH")
#'
#'
#' ### BLAST Reciprocal Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "RBH")
#'
#'
#' ### BLAST Best Hit
#'
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein",
#'           ortho_detection = "BH")
#' }
#' @seealso \code{\link{diamond_rec}}, \code{\link{diamond_best}}, \code{\link{blast_rec}}, \code{\link{dNdS}}
#' @export

orthologs <- function(query_file,
                      subject_files,
                      seq_type           = "protein",
                      outgroup_file      = NULL,
                      eval               = "1E-5",
                      format             = "fasta",
                      ortho_detection    = "DIAMOND_RBH",
                      sensitivity_mode   = "fast",
                      max.target.seqs    = 10000,
                      delete_corrupt_cds = FALSE,
                      cdd.path           = NULL,
                      path               = NULL,
                      add_params         = NULL,
                      diamond_params     = NULL,
                      comp_cores         = 1,
                      quiet              = FALSE,
                      clean_folders      = FALSE) {

        if (!is.element(ortho_detection,
                        c("DIAMOND_RBH", "DIAMOND_BH", "BH", "RBH", "Orthofinder2")))
                stop(
                        paste0(
                                "Unsupported orthology detection method: '", ortho_detection, "'. ",
                                "Choose one of: 'DIAMOND_RBH', 'DIAMOND_BH', 'RBH', 'BH', 'Orthofinder2'."
                        ),
                        call. = FALSE
                )

        i <- NULL

        if (seq_type == "cds") {
                f_sep <- .Platform$file.sep

                filename_qry <-
                        unlist(strsplit(
                                query_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_qry <- filename_qry[length(filename_qry)]

                write.proteome(
                        proteome  = cds2aa(
                                read.cds(query_file, format = format, delete_corrupt_cds = delete_corrupt_cds)
                        ),
                        file.name = file.path(tempdir(), paste0(filename_qry, "_translated.fasta"))
                )

                if (length(subject_files) > 1) {
                        subj_short.names <- vector("character", length(subject_files))

                        for (organism in seq_along(subject_files)) {
                                short.name <-
                                        unlist(
                                                strsplit(
                                                        subject_files[organism],
                                                        f_sep,
                                                        fixed = FALSE,
                                                        perl = TRUE,
                                                        useBytes = FALSE
                                                )
                                        )
                                short.name <- short.name[length(short.name)]
                                subj_short.names[organism] <- short.name

                                write.proteome(
                                        proteome  = cds2aa(
                                                read.cds(subject_files[organism], format = format,
                                                         delete_corrupt_cds = delete_corrupt_cds),
                                                delete_corrupt_cds = delete_corrupt_cds
                                        ),
                                        file.name = file.path(
                                                tempdir(),
                                                paste0(short.name, "_translated.fasta")
                                        )
                                )
                        }

                        subject_files <-
                                file.path(tempdir(),
                                          paste0(subj_short.names, "_translated.fasta"))

                } else {
                        filename_subj <-
                                unlist(
                                        strsplit(
                                                subject_files,
                                                f_sep,
                                                fixed = FALSE,
                                                perl = TRUE,
                                                useBytes = FALSE
                                        )
                                )
                        filename_subj <- filename_subj[length(filename_subj)]

                        write.proteome(
                                proteome  = cds2aa(
                                        read.cds(
                                                subject_files,
                                                format = format,
                                                delete_corrupt_cds = delete_corrupt_cds
                                        ),
                                        delete_corrupt_cds = delete_corrupt_cds
                                ),
                                file.name = file.path(tempdir(),
                                                      paste0(filename_subj, "_translated.fasta"))
                        )

                        subject_files <-
                                file.path(tempdir(),
                                          paste0(filename_subj, "_translated.fasta"))
                }

                query_file <-
                        file.path(tempdir(), paste0(filename_qry, "_translated.fasta"))
        }

        # ---- DIAMOND2 reciprocal best hit (default) ----
        if (ortho_detection == "DIAMOND_RBH") {

                if (length(subject_files) > 1)
                        stop("The DIAMOND2 reciprocal best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        diamond_rec(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                sensitivity_mode   = sensitivity_mode,
                                max.target.seqs    = max.target.seqs,
                                diamond_params     = diamond_params,
                                seq_type           = seq_type,
                                format             = format,
                                clean_folders      = clean_folders
                        )
                )
        }

        # ---- DIAMOND2 best hit ----
        if (ortho_detection == "DIAMOND_BH") {

                if (length(subject_files) > 1)
                        stop("The DIAMOND2 best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        diamond_best(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                sensitivity_mode   = sensitivity_mode,
                                max.target.seqs    = max.target.seqs,
                                diamond_params     = diamond_params,
                                seq_type           = seq_type,
                                format             = format,
                                clean_folders      = clean_folders
                        )
                )
        }

        # ---- BLAST best hit ----
        if (ortho_detection == "BH") {

                if (length(subject_files) > 1)
                        stop("The BLAST best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        blast_best(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                blast_params       = add_params,
                                seq_type           = seq_type,
                                format             = format
                        )
                )

                if (clean_folders)
                        clean_all_folders(file.path(tempdir(), "_blast_db"))
        }

        # ---- BLAST reciprocal best hit ----
        if (ortho_detection == "RBH") {

                if (length(subject_files) > 1)
                        stop("The BLAST reciprocal best hit method is only defined for pairwise comparisons.",
                             call. = FALSE)

                ortho_tbl <- data.table::copy(
                        blast_rec(
                                query_file         = query_file,
                                subject_file       = subject_files,
                                path               = path,
                                delete_corrupt_cds = delete_corrupt_cds,
                                comp_cores         = comp_cores,
                                eval               = eval,
                                blast_params       = add_params,
                                seq_type           = seq_type,
                                format             = format
                        )
                )

                if (clean_folders)
                        clean_all_folders(file.path(tempdir(), "_blast_db"))
        }

        # ---- OrthoFinder2 (under development) ----
        if (ortho_detection == "Orthofinder2") {
                message("The 'Orthofinder2' option is currently under development and will be available soon.")
                return(invisible(NULL))
        }

        return(ortho_tbl)
}
