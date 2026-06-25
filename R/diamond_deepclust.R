#' @title Cluster sequences using \code{diamond deepclust}
#' @description This function clusters a set of protein (or CDS) sequences using the
#' DIAMOND2 \code{deepclust} algorithm. It prepares a DIAMOND2 database via
#' \code{\link{set_diamond}}, runs \code{diamond deepclust}, and returns a
#' two-column tibble mapping each representative sequence to its cluster members.
#' The result is also written to a TSV file.
#' @param input_file a character string specifying the path to the input sequence
#' file (query organism). Supported sequence types are controlled by \code{seq_type}.
#' @param seq_type a character string specifying the sequence type stored in the
#' input file. Options are: \code{"cds"}, \code{"protein"}, or \code{"dna"}.
#' In case of \code{"cds"}, sequences are translated to protein sequences before
#' clustering. Default is \code{seq_type = "protein"}.
#' @param format a character string specifying the file format of the sequence file.
#' Default is \code{format = "fasta"}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with
#' corrupt base triplets should be removed from the input file. Default is
#' \code{delete_corrupt_cds = TRUE}.
#' @param path a character string specifying the path to the DIAMOND2 executable
#' (in case it is not on the system \code{PATH}).
#' @param comp_cores a numeric value specifying the number of CPU threads to use.
#' Passed to \code{--threads}. Default is \code{comp_cores = 1}.
#' @param mutual_cover a numeric value (percentage) specifying the minimum mutual
#' coverage of the cluster member and representative sequence. Passed to
#' \code{--mutual-cover}. Default is \code{NULL} (DIAMOND2 default applies).
#' @param approx_id a numeric value (percentage) specifying the minimum approximate
#' identity to report an alignment or to cluster sequences. Passed to
#' \code{--approx-id}. Default is \code{NULL} (DIAMOND2 default applies).
#' @param eval a numeric value specifying the maximum E-value to report alignments.
#' Passed to \code{--evalue}. Default is \code{NULL} (DIAMOND2 default of 0.001 applies).
#' @param diamond_params a character string of additional DIAMOND2 parameters to
#' append to the \code{deepclust} call. Default is \code{NULL}.
#' @param save.output a character string specifying the path where the output TSV
#' file should be saved. E.g. \code{save.output = getwd()} to save in the current
#' working directory. Default is \code{NULL} (file is only written to a temporary
#' directory).
#' @param quiet a logical value indicating whether DIAMOND2 should run in quiet
#' mode. Default is \code{quiet = TRUE}.
#' @details
#' This function wraps \code{diamond deepclust}, the graph-based sequence clustering
#' algorithm available in DIAMOND2 v2.1.0 and later. The output is a two-column
#' tab-separated file where the first column contains the representative (centroid)
#' sequence accession and the second column contains the cluster member accession.
#' A sequence that is its own representative will appear with itself in both columns.
#'
#' The function uses \code{\link{set_diamond}} internally to handle CDS translation
#' and DIAMOND2 database creation.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments
#' at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki
#' @examples \dontrun{
#' # Cluster CDS sequences (translated to protein internally)
#' diamond_deepclust(
#'   input_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr')
#' )
#'
#' # Cluster protein sequences directly
#' diamond_deepclust(
#'   input_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'   seq_type   = "protein"
#' )
#'
#' # Apply identity and coverage thresholds
#' diamond_deepclust(
#'   input_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   approx_id    = 50,
#'   mutual_cover = 80,
#'   eval         = 1e-5,
#'   comp_cores   = 4
#' )
#'
#' # Save the cluster TSV to the current working directory
#' diamond_deepclust(
#'   input_file  = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'   save.output = getwd()
#' )
#' }
#' @return A \code{\link[tibble]{tibble}} with two columns:
#' \itemize{
#'   \item \code{representative_id} — accession of the cluster representative sequence.
#'   \item \code{member_id} — accession of the cluster member sequence.
#' }
#' @seealso \code{\link{diamond}}, \code{\link{set_diamond}}, \code{\link{diamond_best}}, \code{\link{diamond_rec}}
#' @export
diamond_deepclust <- function(
        input_file,
        seq_type           = "protein",
        format             = "fasta",
        delete_corrupt_cds = TRUE,
        path               = NULL,
        comp_cores         = 1,
        mutual_cover       = NULL,
        approx_id          = NULL,
        eval               = NULL,
        diamond_params     = NULL,
        save.output        = NULL,
        quiet              = TRUE) {

        is_installed_diamond(diamond_exec_path = path)

        if (is.null(path)) {
                message("Running ",
                        system("diamond --version", intern = TRUE)[1],
                        " ...")
        } else {
                message("Running ",
                        system(paste0(
                                'export PATH=$PATH:',
                                path, "' ; diamond --version '"), intern = TRUE)[1],
                        " ...")
        }

        # determine the number of cores on a multicore machine
        cores <- parallel::detectCores()
        if (comp_cores > cores)
                stop("You chose more cores than are available on your machine.",
                     call. = FALSE)

        message("Building DIAMOND2 database for deepclust ...")

        # Translate CDS to protein (if needed) and build a DIAMOND2 database
        db_result <- set_diamond(
                file               = input_file,
                seq_type           = seq_type,
                format             = format,
                makedb             = TRUE,
                delete_corrupt_cds = delete_corrupt_cds,
                path               = path,
                comp_cores         = comp_cores,
                quiet              = quiet
        )

        database <- db_result[[2]]

        # Derive a filename stem for the output TSV
        filename <- unlist(
                strsplit(
                        input_file,
                        .Platform$file.sep,
                        fixed    = FALSE,
                        perl     = TRUE,
                        useBytes = FALSE
                )
        )
        filename <- filename[length(filename)]

        output <- paste0("deepclust_", filename, ".tsv")

        if (!file.exists(file.path(tempdir(), "_blast_db"))) {
                dir.create(file.path(tempdir(), "_blast_db"))
        }

        currwd <- getwd()
        setwd(file.path(tempdir(), "_blast_db"))

        # Build the diamond deepclust command
        deepclust_run <- paste0(
                'diamond deepclust',
                ' --db ', database,
                ' --out ', output,
                ' --threads ', comp_cores
        )

        if (!is.null(approx_id)) {
                deepclust_run <- paste0(
                        deepclust_run,
                        ' --approx-id ', as.numeric(approx_id)
                )
        }

        if (!is.null(mutual_cover)) {
                deepclust_run <- paste0(
                        deepclust_run,
                        ' --mutual-cover ', as.numeric(mutual_cover)
                )
        }

        if (!is.null(eval)) {
                deepclust_run <- paste0(
                        deepclust_run,
                        ' --evalue ', as.numeric(eval)
                )
        }

        if (!is.null(diamond_params)) {
                deepclust_run <- paste0(
                        deepclust_run,
                        ' ', diamond_params
                )
        }

        if (!is.null(path)) {
                deepclust_run <- paste0(
                        'export PATH=$PATH:', path, '; ',
                        deepclust_run
                )
        }

        if (quiet) {
                deepclust_run <- paste0(deepclust_run, ' --quiet')
        }

        message("Running diamond deepclust ...")

        tryCatch({
                system(deepclust_run)
        }, error = function(e) {
                stop(
                        "diamond deepclust did not run correctly.",
                        "\n",
                        "Please check the path to the DIAMOND2 executable and the input file.",
                        "\n",
                        "Error: ", e,
                        call. = FALSE
                )
        })

        tryCatch({
                cluster_table <- data.table::as.data.table(
                        readr::read_tsv(
                                file      = output,
                                col_names = FALSE,
                                col_types = readr::cols(
                                        X1 = readr::col_character(),
                                        X2 = readr::col_character()
                                ),
                                show_col_types = FALSE
                        )
                )

                data.table::setnames(
                        cluster_table,
                        old = c("X1", "X2"),
                        new = c("representative_id", "member_id")
                )

                setwd(file.path(currwd))

                if (!is.null(save.output)) {
                        file.copy(
                                from = file.path(tempdir(), "_blast_db", output),
                                to   = save.output
                        )
                        message("Cluster TSV saved to: ", save.output)
                }

                return(tibble::as_tibble(cluster_table))

        }, error = function(e) {
                stop(
                        "The deepclust output file '", output, "' could not be read correctly.",
                        " Please verify that diamond deepclust produced output.",
                        "\n",
                        "Error: ", e,
                        call. = FALSE
                )
        })
}
