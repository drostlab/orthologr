#' @title Delete the internal folder hierarchy
#' @description This function deletes all internal folders that have been created
#' during pipeline processing. Internally this function uses \code{\link{unlink}}
#' to delete all folders created by the pipline.
#' @param foldernames a character vector storing the folder names that shall be deleted.
#' @author Hajk-Georg Drost
#' @details This function takes a vector storing the names of the folders
#' that shall be deleted, e.g.: \code{clean_all_folders( c("_alignment", "_blast_db", "_calculation") )}.
#' 
#' Since \pkg{orthologr} is a package for pipeline processing based on interface functions,
#' many partial results need to be written to a hard drive to allow subsequent programs
#' to access the output of previous computations. The R core conventions do not favour this
#' behaviour of functions. This would indicate that the \code{clean_folders} argument had to be \code{TRUE}
#' by default in all functions. Nevertheless, this would distract some pipeline functions and therefore
#' \code{clean_folders} = \code{FALSE} by default leaving you with a folder environment that needs
#' to be cleaned by hand, meaning that when using a specific function or pipeline function
#' you must specify \code{clean_folders} = \code{TRUE} in case you want to remove all internal folders
#' after pipeline processing.
#' 
#' 
#' @return This is a void function.
#' @seealso \code{\link{divergence_stratigraphy}}, \code{\link{blast}}, \code{\link{blast_best}}, 
#' \code{\link{blast_rec}}, \code{\link{dNdS}}
#' 
clean_all_folders <- function(foldernames){
        
        if (length(foldernames) > 1) {
                for (i in 1:length(foldernames)) {
                        if (file.exists(foldernames[i])) {
                                unlink(foldernames[i],recursive = TRUE, force = TRUE)
                        }
                }
        } else {
                if (file.exists(foldernames)) {
                        unlink(foldernames,recursive = TRUE, force = TRUE)
                }
        }
}