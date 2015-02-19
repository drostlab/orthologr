#' @title Main Orthology Inference Function
#' @description This function takes nucleotide or protein sequences for a set of organisms 
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' Different orthology inference methods can detect orthologs using multiple subject organisms, e.g. "OrthoMCL", and "PO" (ProteinOrtho).
#' @param outgroup_file a character string specifying the paths to the sequence files of interest (outgroup organisms).
#' This argument is only used by \code{\link{InParanoid}}.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Further methods are: "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), and "OrthoMCL.
#' @param cdd.path path to the cdd database folder (specify when using \code{ortho_detection} = \code{"DELTA"}).
#' @param path a character string specifying the path to the corresponding orthology inference tool.
#' For "BH" and "RBH": path to BLAST, "PO": path to ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL.
#' @param add_params a character string specifying additional parameters that shall be handed to the orthology inference method (tool).
#' Default is \code{add_params} = \code{NULL}.
#' @param detailed_output a boolean value specifying whether a detailed BLAST table shall be returned or only the evalue of the corresponding ortholog pairs. 
#' Default is \code{detailed_output} = \code{TRUE}.
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
#'   \item ProteinOrtho (PO)
#'   \item OrthoMCL (OrthoMCL)
#' }
#' 
#' Alignment based methods:
#' 
#' \itemize{
#'   \item Best Global Alignment (GGSEARCH)
#'   \item Best Local Alignment (SSEARCH)
#' }
#' 
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of orthologous genes
#' in the second column and the amino acid sequences in the third column.
#' @references
#' 
#' BLAST: \url{http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml}
#' 
#' ProteinOrtho: \url{https://www.bioinf.uni-leipzig.de/Software/proteinortho/}
#' 
#' OrthoMCL: \url{http://www.orthomcl.org/orthomcl/}
#' 
#' GGSearch and SSearch: \url{http://fasta.bioch.virginia.edu/fasta_www2/fasta_intro.shtml}
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
#' ### DELTA-BLAST Reciprocal Best Hit
#' 
#' # perform orthology inference using DELTA-BLAST reciprocal best hit
#' # and fasta sequence files storing protein sequences
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "DELTA",
#'           cdd.path        = "path/to/cdd/database/folder")
#'           
#'           
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "DELTA",
#'           cdd.path        = "path/to/cdd/database/folder", 
#'           comp_cores      = 2)          
#'                      
#'           
#' 
#' 
#' ### Orthology Inference using ProteinOrtho
#' 
#' # defining 3 subject organisms: A. lyrata, B. rapa, and T. halophila
#' subject_organisms <- c(system.file('seqs/example_brapa_aa.faa', package = 'orthologr'),
#'                        system.file('seqs/example_alyra_aa.faa', package = 'orthologr'),
#'                        system.file('seqs/example_thalo_aa.faa', package = 'orthologr'))
#' 
#' orthologs(query_file      = system.file('seqs/example_athal_aa.faa', package = 'orthologr'),
#'           subject_files   = subject_organisms,
#'           seq_type        = "protein", 
#'           ortho_detection = "PO")
#'     
#'                 
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = subject_organisms,
#'           seq_type        = "protein", 
#'           ortho_detection = "PO", 
#'           comp_cores      = 2)  
#'           
#'                           
#' ### Orthology Inference using OrthoMCL
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "OrthoMCL")
#'           
#'           
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "OrthoMCL", 
#'           comp_cores      = 2)  
#'
#' 
#' 
#' ### Orthology Inference using GGSearch
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "GGSEARCH")
#'           
#'           
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "GGSEARCH", 
#'           comp_cores      = 2)  
#' 
#'           
#'                               
#' ### Orthology Inference using SSearch
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "SSEARCH")
#'           
#'           
#' # multicore version          
#' orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type        = "protein", 
#'           ortho_detection = "SSEARCH", 
#'           comp_cores      = 2)  
#'
#' }
#' @seealso \code{\link{blast_rec}}, \code{\link{ProteinOrtho}}, \code{\link{OrthoMCL}}, \code{\link{dNdS}}
#' @export
 
orthologs <- function(query_file,
                      subject_files, 
                      seq_type        = "protein",
                      outgroup_file   = NULL, 
                      eval            = "1E-5", 
                      format          = "fasta",
                      ortho_detection = "RBH",
                      cdd.path        = NULL,
                      path            = NULL, 
                      add_params      = NULL,
                      detailed_output = TRUE, 
                      comp_cores      = 1,
                      quiet           = FALSE, 
                      clean_folders   = FALSE){
        
        if(!is.element(ortho_detection, c("BH","RBH","PO","OrthoMCL","GGSEARCH","SSEARCH","DELTA")))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        
        if(ortho_detection == "BH"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        blast_best(query_file      = query_file, 
                                   subject_file    = subject_files, 
                                   path            = path, 
                                   comp_cores      = comp_cores, 
                                   eval            = eval,
                                   blast_params    = add_params, 
                                   seq_type        = seq_type, 
                                   detailed_output = detailed_output, 
                                   format          = format )
                        
                )
                
                
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                
        }
        
        
        if(ortho_detection == "RBH"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        blast_rec( query_file      = query_file, 
                                   subject_file    = subject_files, 
                                   path            = path, 
                                   comp_cores      = comp_cores, 
                                   eval            = eval,
                                   blast_params    = add_params, 
                                   seq_type        = seq_type, 
                                   detailed_output = detailed_output, 
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
                                             subject_file," could not be joined properly to select only the reciprocal best hits.")}
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
        
        
        if(ortho_detection == "OrthoMCL"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        OrthoMCL( query_file      = query_file, 
                                  subject_files   = subject_files, 
                                  orthomcl_path   = path, 
                                  orthomcl_params = add_params,
                                  eval            = eval, 
                                  comp_cores      = comp_cores,
                                  seq_type        = seq_type, 
                                  format          = format )
                        
                        )
                
                
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_OrthoMCL"))
                
        }
        
        
        if(ortho_detection == "GGSEARCH"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        alignmentSearch( query_file      = query_file, 
                                         subject_file   = subject_files,
                                         seq_type        = seq_type,
                                         format          = format,
                                         tool            = "ggsearch",
                                         path            = path, 
                                         comp_cores      = comp_cores,
                                         details         = detailed_output,
                                         clean_folders   = clean_folders,
                                         #eval            = eval 
                                         )
                        
                )
                
        }
        
        
        if(ortho_detection == "SSEARCH"){
                
                
                ortho_tbl <- data.table::copy(
                        
                        alignmentSearch( query_file      = query_file, 
                                         subject_file   = subject_files,
                                         seq_type        = seq_type,
                                         format          = format,
                                         tool            = "ssearch",
                                         path            = path, 
                                         comp_cores      = comp_cores,
                                         details         = detailed_output,
                                         clean_folders   = clean_folders,
                                        #eval            = eval 
                                        )
                        
                )
                
        }
        
#         if(ortho_detection == "IP"){
#                 
#                 
#                 ortho_tbl <- InParanoid(query_file = query_file, subject_file = subject_files, 
#                                    outgroup_file = outgroup_file,ip_path = path,
#                                    seq_type = seq_type, format = format)
#                 
#                 
#                 if(clean_folders)
#                         clean_all_folders("_InParanoid")
#                 
#                 
#         }
#         

        return(ortho_tbl)
        
}