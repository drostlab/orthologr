#' @title Perform a DIAMOND2 search
#' @description This function performs a DIAMOND2 search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param diamond_algorithm a character string specifying the DIAMOND2 algorithm that shall be used, option is currently limited to: \code{diamond_algorithm} = \code{"blastp"}
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' - sensitivity_mode = "faster" : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70
#' - sensitivity_mode = "default" : Default mode. Designed for finding hits of >70
#' - sensitivity_mode = "fast" : fast alignment mode, but least sensitive (default). Designed for finding hits of >70
#' - sensitivity_mode = "mid-sensitive" : fast alignments between the fast mode and the sensitive mode in sensitivity.
#' - sensitivity_mode = "sensitive" : fast alignments, but full sensitivity for hits >40
#' - sensitivity_mode = "more-sensitive" : more sensitive than the sensitive mode.
#' - sensitivity_mode = "very-sensitive" : sensitive alignment mode.
#' - sensitivity_mode = "ultra-sensitive" : most sensitive alignment mode (sensitivity as high as BLASTP).
#' @param eval a numeric value specifying the E-Value cutoff for DIAMOND2 hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param path a character string specifying the path to the DIAMOND2 program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run DIAMOND2 searches.
#' @param diamond_params a character string listing the input parameters that shall be passed to the executing DIAMOND2 program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running DIAMOND2.
#' @param clean_folders a boolean value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the DIAMOND2 output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @param use_blastdb a boolean value specifying whether 
#' Default is \code{use_blastdb} = \code{FALSE}.
#' @details This function provides a fast communication between R and DIAMOND2. It is mainly used as internal functions
#' such as \code{\link{diamond_best}} and \code{\link{diamond_rec}} but can also be used to perform simple DIAMOND2 computations.
#' This function gives the same output as \code{\link{blast}} while being up to 10 000X faster in larger databases.
#'
#' @author Jaruwatana Sodai Lotharukpong
#' @references
#' Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4), 366-368.
#'
#' https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
#' @examples \dontrun{
#' # performing a DIAMOND2 search using diamond blastp (default)
#' diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'
#' # performing a DIAMOND2 search using diamond blastp (default) using amino acid sequences as input file
#' diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'         seq_type     = "protein")
#'
#'
#' # save the DIAMOND2 output table in your current working directory
#' diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'         seq_type     = "protein",
#'         save.output  = getwd())
#'
#' # in case you are working with a multicore machine, you can also run parallel
#' # DIAMOND2 computations using the comp_cores parameter: here with 2 cores
#' diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'         comp_cores   = 2)
#'
#'
#'  # running diamond using additional parameters
#'  diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'          diamond_params = "--max-target-seqs 1")
#'
#'
#' # running diamond using additional parameters and an external diamond path
#'  diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'          diamond_params = "--max-target-seqs 1", path = "path/to/diamond/")
#' }
#'
#' @return A data.table storing the DIAMOND2 hit table returned by DIAMOND2. The format is the same as with BLAST.
#' @seealso \code{\link{diamond_best}}, \code{\link{diamond_rec}}, \code{\link{set_diamond}}, \code{\link{blast}}
#' @export
diamond <- function(
                query_file,
                subject_file,
                seq_type        = "cds",
                format          = "fasta",
                diamond_algorithm = "blastp",
                sensitivity_mode = "fast",
                eval            = "1E-5",
                max.target.seqs = 10000,
                delete_corrupt_cds = TRUE,
                path            = NULL,
                comp_cores      = 1,
                diamond_params    = NULL,
                clean_folders   = FALSE,
                save.output     = NULL,
                database_maker  = "diamond") {
        
        if (!is.element(diamond_algorithm, c("blastp")))
                stop(
                        "Please choose a valid DIAMOND mode. Only 'blastp' is available for this function.",
                        call. = FALSE
                )
        
        if (!is.element(database_maker, c("diamond", "blast")))
                stop("Please choose either: 'diamond' or 'blast' as database_maker.", call. = FALSE)
        
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
        
        message("sensitivity mode: ", sensitivity_mode)
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        aa <- geneids <- NULL
        
        # initialize the DIAMOND search using the previous function for BLAST
        # set_diamond() could be made in the near future. Here is a stopgap.
        if(database_maker == "diamond"){
                message("creating a diamond database")
                query.dt <- set_diamond(
                        file     = query_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds,
                        comp_cores = comp_cores
                )[[1]]
                
                # make a BLASTable/DIAMONDable databse of the subject
                database <- set_diamond(
                        file     = subject_file,
                        seq_type = seq_type,
                        format   = format,
                        makedb   = TRUE ,
                        delete_corrupt_cds = delete_corrupt_cds,
                        comp_cores = comp_cores
                )[[2]]
        }else if(database_maker == "blast"){
                message("creating a blast database")
                query.dt <- set_blast(
                        file     = query_file,
                        seq_type = seq_type,
                        format   = format,
                        delete_corrupt_cds = delete_corrupt_cds
                )[[1]]
                
                # make a BLASTable/DIAMONDable databse of the subject
                database <- set_blast(
                        file     = subject_file,
                        seq_type = seq_type,
                        format   = format,
                        makedb   = TRUE ,
                        delete_corrupt_cds = delete_corrupt_cds
                )[[2]]
        }

        
        filename <-
                unlist(
                        strsplit(
                                query_file,
                                .Platform$file.sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        )
                )
        filename <- filename[length(filename)]
        
        
        # create an internal folder structure for the DIAMOND process
        # we create a BLAST database
        input = paste0("query_", filename, ".fasta")
        # input = "blastinput.fasta"
        output = paste0("blastresult_", filename, ".csv")
        
        
        if (!file.exists(file.path(tempdir(), "_blast_db"))) {
                dir.create(file.path(tempdir(), "_blast_db"))
        }
        
        currwd <- getwd()
        setwd(file.path(tempdir(), "_blast_db"))
        
        # determine the number of cores on a multicore machine
        cores <- parallel::detectCores()
        
        # in case one tries to use more cores than are available
        if (comp_cores > cores)
                stop("You chose more cores than are available on your machine.",
                     call. = FALSE)
        
        #         write_AA <- as.list(query.dt[ ,aa])
        #         name <- query.dt[ ,geneids]
        #
        tryCatch({
                seqinr::write.fasta(
                        sequences = as.list(query.dt[ , aa]),
                        names     = query.dt[ , geneids],
                        nbchar    = 80,
                        open      = "w",
                        file.out  = input
                )
                
        }, error = function(e) {
                stop(
                        "File ",
                        input,
                        " could not be written properly to the internal folder environment.",
                        " Please check the path to ",
                        input,
                        ".",
                        "\n",
                        "Error:",
                        e
                )
        })
        
        # configuring the diamond run in the command line
        # more diamond modes, i.e. blastn, could be added in the future
        # the first step is to make a default run and tag on the
        diamond_run <-  paste0(
                        'diamond blastp --db ',
                        database,
                        ' --query ',
                        input,
                        ' --evalue ',
                        as.numeric(eval),
                        ' --',
                        sensitivity_mode,
                        ' --max-target-seqs ',
                        max.target.seqs,
                        ' --out ',
                        output ,
                        ' --threads ',
                        comp_cores,
                        ' --outfmt 6', ' qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovhsp qcovhsp sstart send slen evalue bitscore score'
                        )
        
        if(!is.null(path)){
                diamond_run <- paste0(
                        'export PATH=$PATH:',
                        path,
                        '; ',
                        diamond_run
                )
        }
        
        if(!is.null(diamond_params)){
                diamond_run <- paste0(
                        diamond_run,
                        ' ',
                        diamond_params
                )
        }

        ## running diamond
        tryCatch({
                system(diamond_run)
        }, error = function(e){
                stop(
                        "Please check the correct path to ",
                        "diamond ",
                        diamond_algorithm,
                        "... the interface call did not work properly.",
                        "\n",
                        "Error:",
                        e
                )
                }
        )
        # additional DIAMOND parameters can be found here:
        # https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
        diamond_table_names <-
                c(
                        "query_id",
                        "subject_id",
                        "perc_identity",
                        "num_ident_matches",
                        "alig_length",
                        "mismatches",
                        "gap_openings",
                        "n_gaps",
                        "pos_match",
                        "ppos",
                        "q_start",
                        "q_end",
                        "q_len",
                        "qcov",
                        "qcovhsp",
                        "s_start",
                        "s_end",
                        "s_len",
                        "evalue",
                        "bit_score",
                        "score_raw"
                )
        # define the colClasses for faster file streaming
        # col_Classes <- c(rep("character", 2), "double", rep("integer", 6), "double", rep("integer", 6), rep("double", 3))
        tryCatch({
                # hit_table <- data.table::fread(
                #         input      = output,
                #         sep        = "\t",
                #         header     = FALSE,
                #         colClasses = col_Classes
                # )
                
                hit_table <-  
                        data.table::as.data.table(
                                readr::read_delim(
                                        file = output, 
                                        delim = "\t", 
                                        col_names = FALSE,
                                        col_types = readr::cols(
                                                "X1" = readr::col_character(),
                                                "X2" = readr::col_character(),
                                                "X3" = readr::col_double(),
                                                "X4" = readr::col_integer(),
                                                "X5" = readr::col_integer(),
                                                "X6" = readr::col_integer(),
                                                "X7" = readr::col_integer(),
                                                "X8" = readr::col_integer(),
                                                "X9" = readr::col_integer(),
                                                "X10" = readr::col_double(),
                                                "X11" = readr::col_integer(),
                                                "X12" = readr::col_integer(),
                                                "X13" = readr::col_integer(),
                                                "X14" = readr::col_double(),
                                                "X15" = readr::col_double(),
                                                "X16" = readr::col_integer(),
                                                "X17" = readr::col_integer(),
                                                "X18" = readr::col_integer(),
                                                "X19" = readr::col_double(),
                                                "X20" = readr::col_number(),   
                                                "X21" = readr::col_double() )))
                
                data.table::setnames(
                        x   = hit_table,
                        old = paste0("X", 1:length(diamond_table_names)),
                        new = diamond_table_names
                )
                
                data.table::setkeyv(hit_table, c("query_id", "subject_id"))
                
                setwd(file.path(currwd))
                
                
                if (clean_folders) {
                        # save the DIAMOND output file to path save.output
                        if (!is.null(save.output))
                                file.copy(file.path(tempdir(), "_blast_db", output),
                                          save.output)
                        
                        clean_all_folders(file.path(tempdir(), "_blast_db"))
                }
                
                if (!clean_folders) {
                        # save the DIAMOND output file to path save.output
                        if (!is.null(save.output))
                                file.copy(file.path(tempdir(), "_blast_db", output),
                                          save.output)
                        
                }
                
                
                hit_table <-
                        tibble::as_tibble(hit_table)
                
                return(hit_table)
                
        }, error = function(e) {
                stop(
                        "File ",
                        output,
                        " could not be read correctly.",
                        " Please check the correct path to ",
                        output,
                        " or whether DIAMOND did write the resulting hit table correctly.",
                        "\n",
                        "Error:",
                        e
                )
        })
}
