#' @title Interface function to InParanoid
#' @description This function takes a query organism and a set of subject organisms 
#' and performs orthology inference using InParanoid as orthology detection program (methodology).
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param outgroup_file a character string specifying the path to the sequence file of interest (outgroup organism).
#' Since the outgroup option in InParanoid is experimental, the default is \code{outgroup_file} = \code{NULL}.
#' @param ip_params a character string specifying additional parameters that shall be handed to the InParanoid call,
#' e.g. \code{ip_params} = \code{""}. See \url{http://inparanoid.sbc.su.se/cgi-bin/faq.cgi}
#' for details.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' There is only the option: "protein". InParanoid can only handle protein sequences stored in fasta files.
#' @param format when using InParanoid you can only specify \code{format} = \code{"fasta"}.
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run InParanoid.
#' @param delete_files a boolean value specifying whether the folder '_InParanoid' that stored the
#' InParanoid output files shall be removed after the analysis. Default is \code{delete_files} = \code{FALSE}.
#' @details This function...
#' @author Hajk-Georg Drost
#' @references
#' 
#' Ostlund G, Schmitt T, Forslund K, Kostler T, Messina DN, Roopra S, Frings O and Sonnhammer ELL.
#' InParanoid 7: new algorithms and tools for eukaryotic orthology analysis. Nucleic Acids Res. 38:D196-D203 (2010)
#'
#' \url{http://inparanoid.sbc.su.se/cgi-bin/index.cgi}
#' 
#' @examples \dontrun{
#' 
#' # finding orthologs between Arabidopsis thaliana and Arabidopsis lyrata genes
#' 
#' }
#' @export
InParanoid <- function(query_file, subject_file, outgroup_file = NULL,ip_params = NULL,eval = "1E-5",
                       seq_type = "protein", format = "fasta", 
                       comp_cores = 1, delete_files = FALSE){
        
        if(format != "fasta")
                stop("InParanoid only supports fasta files.")
        
        if(!is.element(seq_type,c("protein")))
                stop("InParanoid only supports protein sequences. Please choose: seq_type = 'protein'.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!file.exists(paste0("_InParanoid",f_sep))){
                
                dir.create("_InParanoid")
        }
        
        currwd <- getwd()
        setwd(file.path(currwd, "_InParanoid"))
        
        
        
        
        setwd(currwd)
        
        
        
}