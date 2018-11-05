#' @title Main Orthology Inference Function
#' @description This function takes nucleotide or protein sequences for a set of organisms 
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' Different orthology inference methods can detect orthologs using multiple subject organisms, e.g. "OrthoMCL", and "PO" (ProteinOrtho).
#' @param outgroup_file a character string specifying the paths to the sequence files of interest (outgroup organisms).
#' This argument is only used by \code{InParanoid}.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Further methods are: "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), and "OrthoMCL.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet. 
#' @param cdd.path path to the cdd database folder (specify when using \code{ortho_detection} = \code{"DELTA"}).
#' @param path a character string specifying the path to the corresponding orthology inference tool.
#' For "BH" and "RBH": path to BLAST, "PO": path to ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL.
#' @param add_params a character string specifying additional parameters that shall be handed to the orthology inference method (tool).
#' Default is \code{add_params} = \code{NULL}.
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore computations.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @details 
#' This function takes sequence files of a query organism and a subject organism and performs orthology inference
#' using a defined orthology inference method to dectect orthologous genes.
#' 
#' The following interfaces are implemented in the \code{orthologs} function:
#'  
#' BLAST based methods:
#' 
#' \itemize{
#'   \item BLAST best hit (BH) 
#'   \item BLAST reciprocal best hit (RBH)
#'   \item DELTA-BLAST reciprocal best hit (DELTA)
#' }
#' 
#' 
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of orthologous genes
#' in the second column and the amino acid sequences in the third column.
#' @references
#' 
#' BLAST: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
#' 
#' ProteinOrtho: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
#'
#' @examples \dontrun{
#' 
#' 
#' ### BLAST Best Hit
#' 
#' # perform orthology inference using BLAST best hit
#' # and fasta sequence files storing protein sequences
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "BH")
#' 
#' 
#' ### BLAST Reciprocal Best Hit
#' 
#' # perform orthology inference using BLAST reciprocal best hit
#' # and fasta sequence files storing protein sequences
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "RBH")
#'           
#'           
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "RBH", 
#'           comp_cores      = 2)          
#'           
#'           
#'           
#' }
#' @seealso \code{\link{blast_rec}}, \code{\link{dNdS}}
#' @export
 
orthologs <- function(query_file,
                      subject_files, 
                      seq_type        = "protein",
                      outgroup_file   = NULL, 
                      eval            = "1E-5", 
                      format          = "fasta",
                      ortho_detection = "RBH",
                      delete_corrupt_cds = TRUE,
                      cdd.path        = NULL,
                      path            = NULL, 
                      add_params      = NULL,
                      comp_cores      = 1,
                      quiet           = FALSE, 
                      clean_folders   = FALSE){
        
        if(!is.element(ortho_detection,
                       c("BH", "RBH", "PO", "OrthoMCL", "GGSEARCH", "SSEARCH", "DELTA")))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        i <- query_id <- subject_id <- evalue <- NULL
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
                
                write.proteome(proteome  = cds2aa(
                        read.cds(query_file, format = format, delete_corrupt_cds = delete_corrupt_cds)
                ),
                file.name = file.path(tempdir(), paste0(filename_qry, "_translated.fasta")))
                
                if (length(subject_files) > 1) {
                        subj_short.names <- vector("character", length(subject_files))
                        
                        for (organism in 1:length(subject_files)) {
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
                                short.name <-
                                        short.name[length(short.name)]
                                subj_short.names[i] <-  short.name
                                
                                
                                write.proteome(
                                        proteome  = cds2aa(
                                                read.cds(subject_files[organism], format = format, delete_corrupt_cds = delete_corrupt_cds), delete_corrupt_cds = delete_corrupt_cds),
                                        file.name = file.path(
                                                tempdir(),
                                                paste0(short.name, "_translated.fasta")
                                        )
                                )
                        }
                        
                        subject_files <-
                                file.path(tempdir(),
                                          paste0(short.name, "_translated.fasta"))
                        
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
                        filename_subj <-
                                filename_subj[length(filename_subj)]
                        
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
        
        if (ortho_detection == "BH") {
                
                if (length(subject_files) > 1)
                        stop("The BLAST best hit method is only defined for pairwise comparisons.")
                
                ortho_tbl <- data.table::copy(
                        
                        blast_best(query_file      = query_file, 
                                   subject_file    = subject_files, 
                                   path            = path,
                                   delete_corrupt_cds = delete_corrupt_cds,
                                   comp_cores      = comp_cores, 
                                   eval            = eval,
                                   blast_params    = add_params, 
                                   seq_type        = seq_type, 
                                   format          = format )
                        
                )
                if (clean_folders)
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                
        }
        
        
        if(ortho_detection == "RBH"){
                
                if(length(subject_files) > 1)
                        stop("The BLAST best reciprocal hit method is only defined for pairwise comparisons.")
                
                ortho_tbl <- data.table::copy(
                        
                        blast_rec( query_file      = query_file, 
                                   subject_file    = subject_files, 
                                   path            = path,
                                   delete_corrupt_cds = delete_corrupt_cds,
                                   comp_cores      = comp_cores, 
                                   eval            = eval,
                                   blast_params    = add_params, 
                                   seq_type        = seq_type, 
                                   format          = format )
                        
                        )
                
                
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                
        }
        
        
        if(ortho_detection == "DELTA"){
                
                pars <- paste0(add_params," -evalue ",eval," -num_threads ",comp_cores,
                               " -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1")
                
                ortho_tbl_A <- data.table::copy(
                        
                        advanced_blast( query_file      = query_file, 
                                        subject_file    = subject_files,
                                        blast_algorithm = "deltablast",
                                        db_path         = cdd.path,
                                        path            = path, 
                                        blast_params    = pars, 
                                        seq_type        = seq_type,  
                                        format          = format )
                        
                )
                
                ortho_tbl_B <- data.table::copy(
                        
                        advanced_blast( query_file      = subject_files, 
                                        subject_file    = query_file,
                                        blast_algorithm = "deltablast",
                                        db_path         = cdd.path,
                                        path            = path, 
                                        blast_params    = pars, 
                                        seq_type        = seq_type,  
                                        format          = format )
                        
                )
                
                data.table::setnames(ortho_tbl_B, old = c("query_id","subject_id"), new = c("subject_id","query_id"))
                
                tryCatch({       
                        
                        delta_rec_tbl <- dplyr::semi_join(dplyr::tbl_dt(ortho_tbl_A), dplyr::tbl_dt(ortho_tbl_B),
                                                          by = c("query_id","subject_id"))
                        
                        if(detailed_output){
                                
                                return ( delta_rec_tbl )
                        }
                        
                        if(!detailed_output){
                                
                                return ( delta_rec_tbl[ ,list(query_id,subject_id,evalue)] )
                        }
                        
                        
                }, error = function(e){ stop("The BLAST tables resulting from ",query_file, " and ",
                                             subject_files," could not be joined properly to select only the reciprocal best hits.")}
                )
                
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                
        }
        
        if(ortho_detection == "PO"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        ProteinOrtho( query_file    = query_file,
                                      subject_files = subject_files, 
                                      po_path       = path, 
                                      eval          = eval, 
                                      comp_cores    = comp_cores,
                                      po_params     = add_params, 
                                      seq_type      = seq_type, 
                                      format        = format )
                        
                        )
                
               
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_ProteinOrtho"))
                
        }
        
        return(ortho_tbl)
}