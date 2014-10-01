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


set_path <- function(file, add.folder = NULL){
        
        f_sep <- .Platform$file.sep
        file_name <- unlist(strsplit(file,f_sep))
        file_name <- file_name[length(file_name)]
        curr_wd <- unlist(strsplit(getwd(),f_sep))
        wdir <- grepl(" ",curr_wd)
        curr_wd[wdir] <- stringr::str_replace(string = curr_wd[wdir],
                                              replacement = paste0("'",curr_wd[wdir],"'"), 
                                              pattern = curr_wd[wdir])
        
        curr_wd <- paste0(curr_wd,collapse = f_sep)
        
        if(!is.null(add.folder)){
                
                add.folder <- paste0(f_sep,add.folder)
                return(paste0(curr_wd,add.folder,file_name))
                
        } else {
                
                return(paste0(curr_wd,file_name))
        }
        
        
}


