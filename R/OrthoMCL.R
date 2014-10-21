#' @title Interface function to OrthoMCL
#' @description This function takes a query organism and a set of subject organisms 
#' and performs orthology inference using OrthoMCL as orthology detection program (methodology).
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' @param orthomcl_params a character string specifying additional parameters that shall be handed to the OrthoMCL call,
#' e.g. \code{orthomcl_params} = \code{""}. See \url{http://www.orthomcl.org/orthomcl/about.do}
#' for parameter details.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are: "protein", or "dna". OrthoMCL can only handle protein or nucleotide sequences stored in fasta files.
#' @param format when using OrthoMCL you can only specify \code{format} = \code{"fasta"}.
#' @param orthomcl_path a character string specifying the execution path to OrthoMCL.
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run OrthoMCL.
#' @param delete_files a boolean value specifying whether the folder '_OrthoMCL' that stored the
#' OrthoMCL output files shall be removed after the analysis. Default is \code{delete_files} = \code{FALSE}.
#' @details This function...
#' @author Hajk-Georg Drost
#' @references
#' 
#' \url{http://www.orthomcl.org/orthomcl/}
#' @examples \dontrun{
#' 
#' # # finding orthologs between Arabidopsis thaliana and Arabidopsis lyrata genes
#' 
#' }
#' @export

OrthoMCL <- function(query_file, subject_files, orthomcl_params = NULL,eval = "1E-5",
                     seq_type = "protein", format = "fasta", orthomcl_path = NULL, 
                     comp_cores = 1, delete_files = FALSE){
        
        
        
        
        
}