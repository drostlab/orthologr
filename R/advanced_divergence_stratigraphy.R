#' @title A function to perform advanced divergence stratigraphy
#' @description This function uses an advanced methodology to perform divergence stratigraphy.
#' This new methodology needs more computation time, but aims to provide a more
#' accurate and more conservative pipeline for dNdS estimation and divergence stratum assignment.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_set a character string specifying the path to the CDS files of interest (subject organisms).
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform
#'  parallel computations on a multicore machine.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#'  @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out to the console.
#' @param clean_folders a logical value specifying whether the internal folder structure shall be deleted (cleaned) after
#'  processing this function. Default is \code{clean_folders} = \code{FALSE}.
#' @details This function ...
#' @author Hajk-Georg Drost
#' @export

advanced_ds <- function(query_file, subject_set, eval = "1E-5",
                        blast_path = NULL,blast_params = NULL,
                        comp_cores = 1, dnds.threshold = 2, 
                        quiet = FALSE, clean_folders = FALSE ){
        
        
        subj_organisms <- list.files(subject_set)
        
        
        for(i in 1:length(subj_organisms)){
                
                advanced_blast(query_file = query_file,
                               subject_file = subj_organisms[i],
                               seq_type = "cds", blast_algorithm = "blastp",
                               path = blast_path, blast_params = blast_params,
                               makedb_type = "protein")
                                               
        } 
        
        
        
        
        
}