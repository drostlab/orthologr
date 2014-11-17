#' @title Interface function to InParanoid
#' @description This function takes a query organism and a set of subject organisms 
#' and performs orthology inference using InParanoid as orthology detection program (methodology).
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param outgroup_file a character string specifying the path to the sequence file of interest (outgroup organism).
#' Since the outgroup option in InParanoid is experimental, the default is \code{outgroup_file} = \code{NULL}.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' @param ip_path a character string specifying the executable path of InParanoid. Default is \code{ip_path} = \code{NULL}
#' assuming that InParanoid can be called from your default executable PATH.
#' There is only the option: "protein". InParanoid can only handle protein sequences stored in fasta files.
#' @param format when using InParanoid you can only specify \code{format} = \code{"fasta"}.
#' @param delete_files a boolean value specifying whether the folder '_InParanoid' that stored the
#' InParanoid output files shall be removed after the analysis. Default is \code{delete_files} = \code{FALSE}.
#' @details This function aims to provide an interface between \code{R} and the standalone version of the
#' orthology inference program InParanoid.
#' 
#' InParanoid can be downloaded here: \url{http://software.sbc.su.se/cgi-bin/request.cgi?project=inparanoid}.
#' 
#' Note that the standalone version provided by the InParanoid community resembles InParanoid version 4.1 !
#' 
#' To use an outgroup organism (OrthoC) in InParanoid version 4.1, you need to modify the perl script by hand...
#' Within the 'inparanoid.pl' script you need to set the parameter '$use_outgroup = 0;' (default) to '$use_outgroup = 1;',
#' to be able to use the outgroup organism functionality.
#' 
#' In case you changed '$use_outgroup = 1;' in 'inparanoid.pl' you can use the \code{outgroup_file} argument
#' provided by this function (\code{InParanoid}). Otherwise you can only specify \code{query_file} and \code{subject_file} (see examples).
#' 
#' @author Hajk-Georg Drost
#' @references
#' 
#' Ostlund G, Schmitt T, Forslund K, Kostler T, Messina DN, Roopra S, Frings O and Sonnhammer ELL.
#' InParanoid 7: new algorithms and tools for eukaryotic orthology analysis. Nucleic Acids Res. 38:D196-D203 (2010)
#'
#'
#' "InParanoid: A Comprehensive Database of Eukaryotic Orthologs"
#' O'Brien Kevin P, Remm Maido and Sonnhammer Erik L.L
#' Nucleic Acids Res. 33:D476-D480 (2005)
#' 
#' 
#' "Automatic clustering of orthologs and in-paralogs from pairwise species comparisons"
#' Maido Remm, Christian E. V. Storm, and Erik L. L. Sonnhammer
#' J. Mol. Biol. 314:1041-1052 (2001)
#' 
#' 
#' \url{http://inparanoid.sbc.su.se/cgi-bin/index.cgi}
#' 
#' @examples \dontrun{
#' 
#' # finding orthologs between Arabidopsis thaliana and Arabidopsis lyrata genes
#' InParanoid(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'              subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'              seq_type = "protein",format = "fasta")
#'              
#'              
#' }
#' @return a list storing the ortholog groups returned by InParanoid.
#' @seealso \code{\link{orthologs}}, \code{\link{dNdS}}
#' @export
InParanoid <- function(query_file, subject_file, 
                       outgroup_file = NULL,ip_path = NULL,
                       seq_type = "protein", format = "fasta", 
                       delete_files = FALSE){
        
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
        
        
        if(is.null(ip_path)){
                
                if(is.null(outgroup_file)){
                        
                        system(paste0("perl inparanoid.pl ",query_file," ",subject_file))
                                
                } else {
                                
                        system(paste0("perl inparanoid.pl ",query_file," ",subject_file," ",outgroup_file))
                        
                }
                        
                
        } else {
                
                if(is.null(outgroup_file)){
                        
                        system(paste0("perl ",ip_path,f_sep,"inparanoid.pl ",query_file," ",subject_file))
                
                } else {
                        
                        system(paste0("perl ",ip_path,f_sep,"inparanoid.pl ",query_file," ",subject_file," ",outgroup_file))
                        
                }
                                
        }
        
        setwd(currwd)
        
        tryCatch(
                 {
                         # get the query file name
                         query_name <- unlist(strsplit(query_file,f_sep))
                         query_name <- query_name[length(query_name)]
                         
                         # get the subject file name
                         subject_name <- unlist(strsplit(subject_file,f_sep))
                         subject_name <- subject_name[length(subject_name)]
                         
                         if(!is.null(outgroup_file)){
                                 
                                 # get the outgroup file name
                                 outgroup_name <- unlist(strsplit(outgroup_file,f_sep))
                                 outgroup_name <- outgroup_name[length(outgroup_name)]
                                 
                         }
                         
                         
                         if(is.null(outgroup_file)){
                                 
                                 ip_tbl_name <- paste0("table.",query_name,"-",subject_name)
                                 
                         } else {
                                 
                                 ip_tbl_name <- paste0("table.",query_name,"-",subject_name,"-",outgroup_name)
                         }
                         
                         
                         # read InParanoid output
                         InParanoid_list <- sapply(readLines(paste0("_InParanoid",f_sep,ip_tbl_name)),function(x) unlist(strsplit(x,"\t")))
                         
#                          # store the output table of InParanoid
#                          InParanoid_tbl <- data.table::fread(paste0("_InParanoid",f_sep,ip_tbl_name,sep = "\t",header = FALSE, skip = 1)
#                          
#                          if(is.null(outgroup_file)){
#                                  
#                                  data.table::setnames(InParanoid_tbl,old = paste0("V",1:dim(InParanoid_tbl)[2]),
#                                                       new = c("OrthoID","Score","OrthoA","bootstrap_support_A","OrthoB","bootstrap_support_B"))
#                                  
#                                  
#                          } else {
#                                  
#                                  data.table::setnames(InParanoid_tbl,old = paste0("V",1:dim(InParanoid_tbl)[2]),
#                                                       new = c("OrthoID","Score","OrthoA","bootstrap_support_A","OrthoB","bootstrap_support_B","OrthoC","bootstrap_support_C"))
#                                  
#                          }
#                          
#                          data.table::setkeyv(InParanoid_tbl,c("OrthoA","OrthoB"))
#                          
#                          if(delete_files)
#                                  unlink("_InParanoid",recursive = TRUE, force = TRUE)
#                          
#                          
#                          return(InParanoid_tbl)
                  return(InParanoid_list)
        
 
                  }, error = function() stop(paste("The InParanoid_tbl interface call did not terminate properly.",
                                             "Please make sure you passed all parameters correctly to InParanoid_tbl.",sep="\n"))
        )

        
}





