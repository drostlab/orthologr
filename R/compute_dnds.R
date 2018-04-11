#' @title Compute dNdS Values For A Given Pairwise Alignment
#' @description This function takes a vector containing the
#' query amino acid sequence, subject amino acid sequence, query CDS sequence, and subject CDS sequence
#' and then runs the following pipieline:
#' 
#' \itemize{
#' 
#' \item 1) Multiple-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' \item 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' \item 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' }
#' @param complete_tbl a data.table object storing the query_id, subject_id, query_cds (sequence), 
#' subject_cds (sequence), query_aa (sequence), and subject_aa (sequence) of the organisms that shall be compared.
#' @param aa_aln_tool a character string specifying the multiple alignment tool that shall be used for pairwise protein alignments.
#' @param aa_aln_path a character string specifying the path to the corresponding multiple alignment tool.
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type} = "multiple" or \code{aa_aln_type} = \code{"pairwise"}.
#' Default is \code{aa_aln_type} = "multiple".
#' @param aa_aln_params a character string specifying additional parameters that shall be passed to the multiple alignment system call.
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used for codon alignments. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes the amino acid and CDS sequences two orthologous genes
#' and writes the corresponding amino acid and CDS sequences as fasta file into
#' the internal folder environment. The resulting fasta files (two files) store the 
#' amino acid sequence of the query_id and subject_id (file one) and the CDS sequence of
#' the query_id and subject_id (file two). These fasta files are then used to pass through the following pipeline:
#' 
#' 1) Multiple-Alignment or Pairwise-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' @references \url{http://www.r-bloggers.com/the-wonders-of-foreach/}
#' @seealso \code{\link{multi_aln}}, \code{\link{substitutionrate}}, \code{\link{dNdS}}
#' @import foreach
#' @import data.table
compute_dnds <- function(complete_tbl,
                         aa_aln_type     = "multiple", 
                         aa_aln_tool     = "clustalw", 
                         aa_aln_path     = NULL,
                         aa_aln_params   = NULL, 
                         codon_aln_tool  = "pal2nal",
                         dnds_est.method = "YN", 
                         kaks_calc_path  = NULL, 
                         quiet           = FALSE, 
                         comp_cores      = 1, 
                         clean_folders   = FALSE){
        
        if (comp_cores > parallel::detectCores())
                stop("You assigned more cores to the comp_cores argument than are availible on your machine.")
        
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- query_cds <- subject_cds <- query_aa <- subject_aa <- NULL
        
        multicore <- (comp_cores > 1)
        
        orthologs_names <- vector(mode = "list", length = 2)
        cds_seqs <- vector(mode = "list", length = 2)
        aa_seqs <- vector(mode = "list", length = 2)
        cds_session_fasta <- vector(mode = "character", length = 1)
        aa_session_fasta <- vector(mode = "character", length = 1)
        aa_aln_name <- vector(mode = "character", length = 1)
        
        if (!multicore)
                dNdS_values <- vector(mode = "list", length = nrow(complete_tbl))
        
        if (!file.exists(file.path(tempdir(), "_alignment"))) {
                dir.create(file.path(tempdir(), "_alignment"))
        }
        
        if (!file.exists(file.path(tempdir(), "_alignment", "orthologs"))) {
                dir.create(file.path(tempdir(), "_alignment", "orthologs"))
        }
        
        if (!file.exists(file.path(tempdir(), "_alignment", "orthologs", "CDS"))) {
                dir.create(file.path(tempdir(), "_alignment", "orthologs", "CDS"))
        }
        
        if (!file.exists(file.path(tempdir(), "_alignment", "orthologs", "AA"))) {
                dir.create(file.path(tempdir(), "_alignment", "orthologs", "AA"))
        }
        
        
        if (!multicore) {
                for (i in 1:nrow(complete_tbl)) {
                        dNdS_values[i] <- list((function(i) {
                                ### Perform the sampling process in parallel
                                #         dNdS_values <- foreach::foreach(i = 1:nrow(complete_tbl),.combine="rbind") %dopar%{
                                
                                # storing the query gene id and subject gene id of the orthologous gene pair
                                orthologs_names <-
                                        list(complete_tbl[i , query_id], complete_tbl[i, subject_id])
                                
                                # storing the query CDS sequence and subject CDS sequence of the orthologous gene pair
                                cds_seqs <-
                                        list(complete_tbl[i, query_cds], complete_tbl[i, subject_cds])
                                
                                # storing the query amino acid sequence and subject amino acid sequence of the orthologous gene pair
                                aa_seqs <-
                                        list(complete_tbl[i, query_aa], complete_tbl[i, subject_aa])
                                
                                
                                #aa_aln_name <- paste0("query_",i,"___","subject_",i)
                                aa_aln_name <- paste0("q", i)
                                
                                # create cds fasta of orthologous gene pair having session name: 'aa_aln_name'
                                cds_session_fasta <-
                                        file.path(
                                                tempdir(),
                                                "_alignment",
                                                "orthologs",
                                                "CDS",
                                                paste0(aa_aln_name, "_cds.fasta")
                                        )
                                
                                seqinr::write.fasta(
                                        sequences = cds_seqs,
                                        names     = orthologs_names,
                                        file.out  = cds_session_fasta
                                )
                                
                                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                                aa_session_fasta <-
                                        file.path(
                                                tempdir(),
                                                "_alignment",
                                                "orthologs",
                                                "AA",
                                                paste0(aa_aln_name, "_aa.fasta")
                                        )
                                
                                seqinr::write.fasta(
                                        sequences = aa_seqs,
                                        names     = orthologs_names,
                                        file.out  = aa_session_fasta
                                )
                                
                                # which multi_aln tool should get the parameters
                                #multi_aln_tool_params <- paste0(aa_aln_tool,".",params)
                                
                                
                                #                 pairwise_aln <- Biostrings::pairwiseAlignment(aa_seqs[[1]],aa_seqs[[2]], type = "global")
                                #
                                #                 Biostrings::writePairwiseAlignments(pairwise_aln, block.width = 60)
                                #
                                
                                if (aa_aln_type == "multiple") {
                                        # align aa -> <aa_aln_tool>.aln
                                        multi_aln(
                                                file           = aa_session_fasta,
                                                tool           = aa_aln_tool,
                                                get_aln        = FALSE,
                                                multi_aln_name = aa_aln_name,
                                                path           = aa_aln_path,
                                                quiet          = quiet
                                        )
                                        
                                        aa_aln_session_name <-
                                                paste0(aa_aln_name,
                                                       "_",
                                                       aa_aln_tool,
                                                       ".aln")
                                        
                                        # align codon -> cds.aln
                                        codon_aln(
                                                file_aln       = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "multi_aln",
                                                        paste0(aa_aln_session_name)
                                                ),
                                                file_nuc       = cds_session_fasta,
                                                tool           = codon_aln_tool,
                                                format         = "fasta",
                                                codon_aln_name = aa_aln_name,
                                                get_aln        = FALSE,
                                                quiet          = quiet
                                        )
                                }
                                
                                if (aa_aln_type == "pairwise") {
                                        # align aa -> <aa_aln_tool>.aln
                                        pairwise_aln(
                                                file              = aa_session_fasta,
                                                tool              = aa_aln_tool,
                                                seq_type          = "protein",
                                                get_aln           = FALSE,
                                                pairwise_aln_name = aa_aln_name,
                                                path              = aa_aln_path,
                                                quiet             = quiet
                                        )
                                        
                                        aa_aln_session_name <-
                                                paste0(aa_aln_name,
                                                       "_",
                                                       aa_aln_tool,
                                                       "_AA.aln")
                                        
                                        # align codon -> cds.aln
                                        codon_aln(
                                                file_aln       = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "pairwise_aln",
                                                        paste0(aa_aln_session_name)
                                                ),
                                                file_nuc       = cds_session_fasta,
                                                tool           = codon_aln_tool,
                                                format         = "fasta",
                                                codon_aln_name = aa_aln_name,
                                                get_aln        = FALSE,
                                                quiet          = quiet
                                        )
                                }
                                
                                codon_aln_session_name <-
                                        paste0(aa_aln_name,
                                               "_",
                                               codon_aln_tool,
                                               ".aln")
                                
                                # compute kaks
                                dNdS.table <-
                                        substitutionrate(
                                                file           = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "codon_aln",
                                                        codon_aln_session_name
                                                ),
                                                subst_name     = aa_aln_name,
                                                kaks_calc_path = kaks_calc_path,
                                                est.method     = dnds_est.method,
                                                quiet          = quiet
                                        )
                                
                                # res <- dplyr::inner_join(dtplyr::tbl_dt(dNdS.table), dtplyr::tbl_dt(complete_tbl), by = "query_id")
                                # return(res)
                                return(dNdS.table)
                        })(i))
                }
        }
        
        if (multicore) {
                ### Parallellizing the sampling process using the 'doParallel' and 'parallel' package
                ### register all given cores for parallelization
                clust <- parallel::makeCluster(comp_cores)
                
                ### Perform the sampling process in parallel
                dNdS_values <-
                        parallel::parLapply(clust, seq_len(nrow(complete_tbl)), function(i) {
                                
                                # storing the query gene id and subject gene id of the orthologous gene pair
                                orthologs_names <-
                                        list(complete_tbl[i , query_id], complete_tbl[i, subject_id])
                                
                                # storing the query CDS sequence and subject CDS sequence of the orthologous gene pair
                                cds_seqs <-
                                        list(complete_tbl[i, query_cds], complete_tbl[i, subject_cds])
                                
                                # storing the query amino acid sequence and subject amino acid sequence of the orthologous gene pair
                                aa_seqs <-
                                        list(complete_tbl[i, query_aa], complete_tbl[i, subject_aa])
                                
                                #aa_aln_name <- paste0("query_",i,"___","subject_",i)
                                aa_aln_name <- paste0("q", i)
                                
                                # create cds fasta of orthologous gene pair having session name: 'aa_aln_name'
                                cds_session_fasta <-
                                        file.path(
                                                tempdir(),
                                                "_alignment",
                                                "orthologs",
                                                "CDS",
                                                paste0(aa_aln_name, "_cds.fasta")
                                        )
                                
                                seqinr::write.fasta(
                                        sequences = cds_seqs,
                                        names     = orthologs_names,
                                        file.out  = cds_session_fasta
                                )
                                
                                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                                aa_session_fasta <-
                                        file.path(
                                                tempdir(),
                                                "_alignment",
                                                "orthologs",
                                                "AA",
                                                paste0(aa_aln_name, "_aa.fasta")
                                        )
                                
                                seqinr::write.fasta(
                                        sequences = aa_seqs,
                                        names     = orthologs_names,
                                        file.out  = aa_session_fasta
                                )
                                
                                # which multi_aln tool should get the parameters
                                #multi_aln_tool_params <- paste0(aa_aln_tool,".",params)
                                
                                #                 pairwise_aln <- Biostrings::pairwiseAlignment(aa_seqs[[1]],aa_seqs[[2]], type = "global")
                                #
                                #                 Biostrings::writePairwiseAlignments(pairwise_aln, block.width = 60)
                                #
                                
                                if (aa_aln_type == "multiple") {
                                        # align aa -> <aa_aln_tool>.aln
                                        multi_aln(
                                                file           = aa_session_fasta,
                                                tool           = aa_aln_tool,
                                                get_aln        = FALSE,
                                                multi_aln_name = aa_aln_name,
                                                path           = aa_aln_path,
                                                quiet          = quiet
                                        )
                                        
                                        aa_aln_session_name <-
                                                paste0(aa_aln_name, "_", aa_aln_tool, ".aln")
                                        
                                        # align codon -> cds.aln
                                        codon_aln(
                                                file_aln       = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "multi_aln",
                                                        aa_aln_session_name
                                                ),
                                                file_nuc       = cds_session_fasta,
                                                tool           = codon_aln_tool,
                                                format         = "fasta",
                                                codon_aln_name = aa_aln_name,
                                                get_aln        = FALSE,
                                                quiet          = quiet
                                        )
                                }
                                
                                if (aa_aln_type == "pairwise") {
                                        # align aa -> <aa_aln_tool>.aln
                                        pairwise_aln(
                                                file              = aa_session_fasta,
                                                tool              = aa_aln_tool,
                                                seq_type          = "protein",
                                                get_aln           = FALSE,
                                                pairwise_aln_name = aa_aln_name,
                                                path              = aa_aln_path,
                                                quiet             = quiet
                                        )
                                        
                                        aa_aln_session_name <-
                                                paste0(aa_aln_name,
                                                       "_",
                                                       aa_aln_tool,
                                                       "_AA.aln")
                                        
                                        # align codon -> cds.aln
                                        codon_aln(
                                                file_aln       = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "pairwise_aln",
                                                        aa_aln_session_name
                                                ),
                                                file_nuc       = cds_session_fasta,
                                                tool           = codon_aln_tool,
                                                format         = "fasta",
                                                codon_aln_name = aa_aln_name,
                                                get_aln        = FALSE,
                                                quiet          = quiet
                                        )
                                }
                                
                                codon_aln_session_name <-
                                        paste0(aa_aln_name, "_", codon_aln_tool, ".aln")
                                
                                # compute kaks
                                dNdS.table <-
                                        substitutionrate(
                                                file           = file.path(
                                                        tempdir(),
                                                        "_alignment",
                                                        "codon_aln",
                                                        codon_aln_session_name
                                                ),
                                                subst_name     = aa_aln_name,
                                                kaks_calc_path = kaks_calc_path,
                                                est.method     = dnds_est.method,
                                                quiet          = quiet
                                        )
                                
                                # res <- dplyr::inner_join(dtplyr::tbl_dt(dNdS.table), dtplyr::tbl_dt(complete_tbl), by = "query_id")
                                # return(res)
                                return(dNdS.table)

                        })
                
                parallel::stopCluster(clust)
        }
        
        
        dNdS_tbl <- data.table::as.data.table(do.call(rbind, dNdS_values))
        
        setkeyv(dNdS_tbl, c("query_id", "subject_id"))
        
        
        if (clean_folders)
                clean_all_folders(c(
                        file.path(tempdir(), "_alignment"),
                        file.path(tempdir(), "_blast_db"),
                        file.path(tempdir(), "_calculation")
                ))
        
        # returning the dNdS table as data.table object
        return(dNdS_tbl)
        
}
