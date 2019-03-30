#' @title Perform a BLAST+ search
#' @description This function performs a BLAST+ search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, option is: \code{blast_algorithm} = \code{"blastp"}
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param remote a boolean value specifying whether a remote BLAST search shall be performed.
#' In case \code{remote} = \code{TRUE}, please specify the \code{db} argument. This feature is very experimental,
#' since a query of only a few genes against NCBI nr database, can consume a lot of time and might cause
#' response delay crashes in R.
#' @param db a character string specifying the NCBI data base that shall be queried using remote BLAST.
#' This parameter must be specified when \code{remote} = \code{TRUE} and is \code{NULL} by default.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param clean_folders a boolean value specifying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the BLAST output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @details This function provides a fast communication between R and BLAST+. It is mainly used as internal functions
#' such as \code{\link{blast_best}} and \code{\link{blast_rec}} but can also be used to perform simple BLAST computations.
#'
#' Note, that this function isn't as flexible as \code{\link{advanced_blast}}.
#'
#' When using \code{remote} = \code{TRUE}, make sure you specify the \code{db} argument.
#' The following databases can be chosen:
#'
#' \code{db}
#' \itemize{
#' \item "nr"
#' \item "plaza"
#' }
#'
#' Note: When working with remote BLAST, make sure you don't submit too large jobs due to the
#' BLAST query conventions! All in all this functionality is still very experimental and can cause
#' problems due to time out errors when submitting too large queries!
#'
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @references
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#'
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schaeffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
#'
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi
#' @examples \dontrun{
#' # performing a BLAST search using blastp (default)
#' blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'
#' # performing a BLAST search using blastp (default) using amino acid sequences as input file
#' blast(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type     = "protein")
#'
#'
#' # save the BLAST output table in your current working directory
#' blast(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type     = "protein",
#'       save.output  = getwd())
#'
#' # in case you are working with a multicore machine, you can also run parallel
#' # BLAST computations using the comp_cores parameter: here with 2 cores
#' blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       comp_cores   = 2)
#'
#'
#'  # running blastp using additional parameters
#'  blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'        blast_params = "-max_target_seqs 1")
#'
#'
#' # running blastp using additional parameters and an external blastp path
#'  blast(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'        blast_params = "-max_target_seqs 1", path = "path/to/blastp/")
#' }
#'
#' @return A data.table storing the BLAST hit table returned by BLAST.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
blast <- function(query_file,
                  subject_file,
                  seq_type        = "cds",
                  format          = "fasta",
                  blast_algorithm = "blastp",
                  eval            = "1E-5",
                  max.target.seqs = 10000,
                  delete_corrupt_cds = TRUE,
                  remote          = FALSE,
                  db              = NULL,
                  path            = NULL,
                  comp_cores      = 1,
                  blast_params    = NULL,
                  clean_folders   = FALSE,
                  save.output     = NULL) {
        if (!is.element(blast_algorithm, c("blastp")))
                stop(
                        "Please choose a valid BLAST mode. Only 'blastp' is available for this function.",
                        call. = FALSE
                )
        
        if (remote & is.null(db))
                stop(
                        "To use the remote option of blast() please specify the 'db' argument, e.g. db = 'nr'",
                        call. = FALSE
                )
        
        if (!is.null(db)) {
                if (!is.element(db, c("nr", "plaza")))
                        stop("Please choose a database that is supported by remote BLAST.",
                             call. = FALSE)
        }
        
        is_installed_blast()
        message("Running ",
                system("blastp -version", intern = TRUE)[1],
                " ...")
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        aa <- geneids <- NULL
        
        # initialize the BLAST search
        query.dt <- set_blast(
                file     = query_file,
                seq_type = seq_type,
                format   = format,
                delete_corrupt_cds = delete_corrupt_cds
        )[[1]]
        
        # make a BLASTable databse of the subject
        database <- set_blast(
                file     = subject_file,
                seq_type = seq_type,
                format   = format,
                makedb   = TRUE ,
                delete_corrupt_cds = delete_corrupt_cds
        )[[2]]
        
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
        
        
        # create an internal folder structure for the BLAST process
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
                        sequences = as.list(query.dt[, aa]),
                        names     = query.dt[, geneids],
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
                        "."
                )
        })
        
        # test whether the connection to BLAST+ works
        tryCatch({
                if (remote) {
                        # use the default parameters when running blastp
                        system(
                                paste0(
                                        "blastp -db ",
                                        db,
                                        " -query ",
                                        input,
                                        " -remote -evalue ",
                                        eval,
                                        " -out ",
                                        output ,
                                        " -outfmt 6"
                                )
                        )
                        
                } else {
                        if (is.null(path)) {
                                if (blast_algorithm == "blastp") {
                                        if (is.null(blast_params)) {
                                                # use the default parameters when running blastp
                                                system(
                                                        paste0(
                                                                "blastp -db ",
                                                                database,
                                                                " -query ",
                                                                input,
                                                                " -evalue ",
                                                                eval,
                                                                " -max_target_seqs ",
                                                                max.target.seqs,
                                                                " -out ",
                                                                output ,
                                                                " -outfmt 6",
                                                                " -num_threads ",
                                                                comp_cores
                                                        )
                                                )
                                        } else {
                                                # add additional parameters when running blastp
                                                system(
                                                        paste0(
                                                                "blastp -db ",
                                                                database,
                                                                " -query ",
                                                                input,
                                                                " -evalue ",
                                                                eval,
                                                                " -max_target_seqs ",
                                                                max.target.seqs,
                                                                " -out ",
                                                                output ,
                                                                " -outfmt 6",
                                                                " -num_threads ",
                                                                comp_cores,
                                                                " ",
                                                                blast_params
                                                        )
                                                )
                                                
                                        }
                                }
                                
                        } else {
                                if (blast_algorithm == "blastp") {
                                        if (is.null(blast_params)) {
                                                # use the default parameters when running blastp
                                                system(
                                                        paste0(
                                                                "export PATH=$PATH:",
                                                                path,
                                                                "; blastp -db ",
                                                                database,
                                                                " -query ",
                                                                input,
                                                                " -evalue ",
                                                                eval,
                                                                " -max_target_seqs ",
                                                                max.target.seqs,
                                                                " -out ",
                                                                output ,
                                                                " -outfmt 6",
                                                                " -num_threads ",
                                                                comp_cores
                                                        )
                                                )
                                        } else {
                                                # add additional parameters when running blastp
                                                system(
                                                        paste0(
                                                                "export PATH=$PATH:",
                                                                path,
                                                                "; blastp -db ",
                                                                database,
                                                                " -query ",
                                                                input,
                                                                " -evalue ",
                                                                eval,
                                                                " -max_target_seqs ",
                                                                max.target.seqs,
                                                                " -out ",
                                                                output ,
                                                                " -outfmt 6",
                                                                " -num_threads ",
                                                                comp_cores,
                                                                " ",
                                                                blast_params
                                                        )
                                                )
                                                
                                                
                                        }
                                }
                        }
                }
                
        }, error = function(e) {
                stop(
                        "Please check the correct path to ",
                        blast_algorithm,
                        "... the interface call did not work properly."
                )
        })
        
        # additional blast parameters can be found here:
        # http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
        blast_table_names <-
                c(
                        "query_id",
                        "subject_id",
                        "perc_identity",
                        "alig_length",
                        "mismatches",
                        "gap_openings",
                        "q_start",
                        "q_end",
                        "s_start",
                        "s_end",
                        "evalue",
                        "bit_score"
                )
        
        # define the colClasses for faster file streaming
        col_Classes <- c(rep("character", 2), rep("numeric", 10))
        
        tryCatch({
                hit_table <- data.table::fread(
                        input      = output,
                        sep        = "\t",
                        header     = FALSE,
                        colClasses = col_Classes
                )
                
                data.table::setnames(
                        x   = hit_table,
                        old = paste0("V", 1:length(blast_table_names)),
                        new = blast_table_names
                )
                
                data.table::setkeyv(hit_table, c("query_id", "subject_id"))
                
                setwd(file.path(currwd))
                
                
                if (clean_folders) {
                        # save the BLAST output file to path save.output
                        if (!is.null(save.output))
                                file.copy(file.path(tempdir(), "_blast_db", output),
                                          save.output)
                        
                        clean_all_folders(file.path(tempdir(), "_blast_db"))
                }
                
                if (!clean_folders) {
                        # save the BLAST output file to path save.output
                        if (!is.null(save.output))
                                file.copy(file.path(tempdir(), "_blast_db", output),
                                          save.output)
                        
                }
                
                
                hit_table <-
                        tibble::as_tibble(dtplyr::tbl_dt(hit_table))
                
                return(hit_table)
                
        }, error = function(e) {
                stop(
                        "File ",
                        output,
                        " could not be read correctly.",
                        " Please check the correct path to ",
                        output,
                        " or whether BLAST did write the resulting hit table correctly."
                )
        })
        
}
