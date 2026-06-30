#' @title Realign sequences within \code{deepclust} clusters
#' @description Runs \code{diamond realign} on an existing set of
#' \code{\link{deepclust}} clusters to obtain per-pair alignment statistics
#' (approximate identity, e-value, bitscore, and query/subject coverages) for
#' every \code{(member, representative)} pair inside each cluster. The result
#' can be passed directly to \code{\link{deepclust_annotate}} via its
#' \code{realign} argument to enrich a cluster profile without re-running
#' \code{diamond deepclust}.
#' @param input_file a character string, or a character vector of file paths,
#' specifying the input sequence file(s) that were used for the original
#' \code{\link{deepclust}} run. When multiple paths are provided they are merged
#' into a single temporary FASTA before the DIAMOND2 database is built (same
#' behaviour as \code{\link{deepclust}}). Both plain and gzip-compressed
#' (\code{.gz}) files are supported.
#' @param clusters either
#' \itemize{
#'   \item a \code{character} string — path to a tab-separated
#'     \code{deepclust} output file whose first column is
#'     \code{representative_id} and second column is \code{member_id}, or
#'   \item a \code{data.frame} / \code{tibble} with columns
#'     \code{representative_id} and \code{member_id} (i.e. the direct return
#'     value of \code{\link{deepclust}} or \code{\link{deepclust_annotate}}).
#' }
#' @param seq_type a character string specifying the sequence type stored in the
#' input file(s). Options are: \code{"cds"}, \code{"protein"}, or \code{"dna"}.
#' In case of \code{"cds"}, sequences are translated to protein sequences before
#' building the database. Default is \code{seq_type = "protein"}.
#' @param format a character string specifying the file format of the sequence
#' file. Default is \code{format = "fasta"}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with
#' corrupt base triplets should be removed from the input file. Only relevant
#' when \code{seq_type = "cds"}. Default is \code{delete_corrupt_cds = TRUE}.
#' @param path a character string specifying the path to the DIAMOND2 executable
#' (in case it is not on the system \code{PATH}).
#' @param comp_cores a numeric value specifying the number of CPU threads to
#' use. Passed to \code{--threads}. Default is \code{comp_cores = 1}.
#' @param realign_params a character string of additional \code{diamond realign}
#' parameters to append to the command. Default is \code{NULL}.
#' @param save.output a character string specifying the path where the output
#' TSV file should be saved. E.g. \code{save.output = getwd()} to save in the
#' current working directory. Default is \code{NULL} (file is only written to a
#' temporary directory).
#' @param quiet a logical value indicating whether DIAMOND2 should run in quiet
#' mode. Default is \code{quiet = TRUE}.
#' @details
#' This function wraps \code{diamond realign}, which re-aligns every cluster
#' member against its assigned representative using the cluster mapping produced
#' by \code{\link{deepclust}}. The output format is BLAST tabular
#' (\code{--outfmt 6}) with the fields:
#' \code{qseqid sseqid approx_pident evalue bitscore qcovhsp scovhsp}.
#'
#' The \code{--clusters} argument consumed by \code{diamond realign} expects a
#' two-column tab-separated file (\code{representative_id TAB member_id}).
#' When a tibble is supplied to \code{clusters}, it is written to a temporary
#' file automatically.
#'
#' Columns in the returned tibble are named to match those used by
#' \code{\link{deepclust}} (\code{representative_id}, \code{member_id}) so that
#' the result joins cleanly onto \code{\link{deepclust_annotate}} output.
#' Note that \code{diamond realign} uses the \emph{representative} as the
#' alignment query (\code{qseqid}) and the \emph{member} as the subject
#' (\code{sseqid}); the returned tibble is named accordingly.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments
#' at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki
#' @examples \dontrun{
#' # 1. Run deepclust first
#' clusters <- deepclust(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   )
#' )
#'
#' # 2. Realign each representative against its cluster members to get alignment stats
#' realign_result <- deepclust_realign(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   clusters = clusters
#' )
#'
#' # 3. Alternatively, supply the path to the deepclust TSV directly
#' realign_result <- deepclust_realign(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   clusters = "/path/to/deepclust_output.tsv"
#' )
#'
#' # 4. Incorporate into deepclust_annotate without re-running deepclust
#' profile <- deepclust_annotate(
#'   input_file = c(
#'     system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'     system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
#'   ),
#'   realign = realign_result
#' )
#' }
#' @return A \code{\link[tibble]{tibble}} with seven columns:
#' \itemize{
#'   \item \code{representative_id} — accession of the cluster representative
#'     (query in the realignment).
#'   \item \code{member_id} — accession of the cluster member (subject in the
#'     realignment).
#'   \item \code{approx_pident} — approximate percentage of identical matches.
#'   \item \code{evalue} — expect value.
#'   \item \code{bitscore} — bit score.
#'   \item \code{qcovhsp} — query coverage per HSP.
#'   \item \code{scovhsp} — subject (representative) coverage per HSP.
#' }
#' @seealso \code{\link{deepclust}}, \code{\link{deepclust_annotate}},
#'   \code{\link{diamond}}, \code{\link{set_diamond}}
#' @export
deepclust_realign <- function(
        input_file,
        clusters,
        seq_type           = "protein",
        format             = "fasta",
        delete_corrupt_cds = TRUE,
        path               = NULL,
        comp_cores         = 1,
        realign_params     = NULL,
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

        # validate that all supplied sequence files exist before doing any work
        missing_files <- input_file[!file.exists(input_file)]
        if (length(missing_files) > 0)
                stop(
                        "The following input file(s) were not found:\n",
                        paste(missing_files, collapse = "\n"),
                        call. = FALSE
                )

        # ensure the shared temp directory exists early
        if (!file.exists(file.path(tempdir(), "_blast_db")))
                dir.create(file.path(tempdir(), "_blast_db"))

        # resolve the clusters argument to a file path
        if (is.character(clusters) && length(clusters) == 1 && file.exists(clusters)) {
                clusters_file <- clusters
        } else if (is.data.frame(clusters)) {
                required_cols <- c("representative_id", "member_id")
                if (!all(required_cols %in% names(clusters)))
                        stop(
                                "When 'clusters' is a data frame it must have columns: ",
                                paste(required_cols, collapse = ", "),
                                call. = FALSE
                        )
                clusters_file <- file.path(tempdir(), "_blast_db",
                                           "deepclust_realign_input.tsv")
                readr::write_tsv(
                        clusters[, required_cols],
                        file      = clusters_file,
                        col_names = FALSE
                )
        } else {
                stop(
                        "'clusters' must be either a path to a deepclust TSV file or a ",
                        "data frame with columns 'representative_id' and 'member_id'.",
                        call. = FALSE
                )
        }

        # merge multiple input files into a single temporary FASTA before
        # building the database (mirrors deepclust behaviour)
        if (length(input_file) > 1) {

                n_files       <- length(input_file)
                seqtype_fasta <- if (seq_type == "protein") "AA" else "DNA"

                message("Merging ", n_files, " input files ...")

                all_seqs  <- list()
                all_names <- character(0)

                for (f in input_file) {
                        seqs <- seqinr::read.fasta(
                                file            = f,
                                seqtype         = seqtype_fasta,
                                as.string       = TRUE,
                                forceDNAtolower = FALSE
                        )
                        all_seqs  <- c(all_seqs,  seqs)
                        all_names <- c(all_names, names(seqs))
                }

                merged_filename <- paste0("deepclust_merged_", n_files, "_files.fasta")
                merged_file     <- file.path(tempdir(), "_blast_db", merged_filename)

                seqinr::write.fasta(
                        sequences = all_seqs,
                        names     = all_names,
                        file.out  = merged_file,
                        nbchar    = 80
                )

                message("Merged ", length(all_names), " sequences into: ", merged_filename)

                input_file <- merged_file
        }

        # derive output filename from the (possibly merged) input file
        filename <- basename(input_file)
        output   <- paste0("deepclust_realign_", filename, ".tsv")

        message("Building DIAMOND2 database for realign ...")

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

        currwd <- getwd()
        setwd(file.path(tempdir(), "_blast_db"))
        on.exit(setwd(currwd), add = TRUE)

        realign_run <- paste0(
                'diamond realign',
                ' --db ',       database,
                ' --clusters ', clusters_file,
                ' --out ',      output,
                ' --threads ',  comp_cores,
                ' --outfmt 6 qseqid sseqid approx_pident evalue bitscore qcovhsp scovhsp'
        )

        if (!is.null(realign_params))
                realign_run <- paste0(realign_run, ' ', realign_params)

        if (!is.null(path))
                realign_run <- paste0('export PATH=$PATH:', path, '; ', realign_run)

        if (quiet)
                realign_run <- paste0(realign_run, ' --quiet')

        message("Running diamond realign ...")

        status <- system(realign_run)
        if (!identical(status, 0L)) {
                stop(
                        "diamond realign exited with non-zero status: ", status,
                        "\nPlease check the path to the DIAMOND2 executable, the input_file, and clusters.",
                        call. = FALSE
                )
        }

        tryCatch({
                recluster_table <- data.table::as.data.table(
                        readr::read_tsv(
                                file      = output,
                                col_names = FALSE,
                                col_types = readr::cols(
                                        X1 = readr::col_character(),
                                        X2 = readr::col_character(),
                                        X3 = readr::col_double(),
                                        X4 = readr::col_double(),
                                        X5 = readr::col_double(),
                                        X6 = readr::col_double(),
                                        X7 = readr::col_double()
                                ),
                                show_col_types = FALSE
                        )
                )

                data.table::setnames(
                        recluster_table,
                        old = paste0("X", 1:7),
                        new = c("representative_id", "member_id",
                                "approx_pident", "evalue", "bitscore",
                                "qcovhsp", "scovhsp")
                )

                setwd(file.path(currwd))

                if (!is.null(save.output)) {
                        file.copy(
                                from = file.path(tempdir(), "_blast_db", output),
                                to   = save.output
                        )
                        message("Recluster TSV saved to: ", save.output)
                }

                return(tibble::as_tibble(recluster_table))

        }, error = function(e) {
                stop(
                        "The deepclust_realign output file '", output,
                        "' could not be read correctly.",
                        " Please verify that diamond realign produced output.",
                        "\n",
                        "Error: ", e,
                        call. = FALSE
                )
        })
}
