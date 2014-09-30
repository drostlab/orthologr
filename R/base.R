#' @title Function to delete the internal folder hierarchy
#' @description This function deletes all internal folders that have been created
#' during pipeline processing. Internally this function uses \code{\link{unlink}}
#' to delete all folders created by the pipline.
#' @author Hajk-Georg Drost
#' @details Following folders are being removed: "_alignment", "_blast_db", and "_calculation".
#' @return Nothing is being returned.
#' @examples \dontrun{
#' 
#' # in case internal folders exist, they are being removed
#' clean_all_folders() 
#' 
#' }
#' @export
clean_all_folders <- function(){
        
        if(file.exists("_alignment")){
                unlink("_alignment",recursive = TRUE, force = TRUE)
        }
        
        if(file.exists("_blast_db")){
                unlink("_blast_db",recursive = TRUE, force = TRUE)
        }
                
        
        if(file.exists("_calculation")){
                unlink("_calculation",recursive = TRUE, force = TRUE)
        }
                
}