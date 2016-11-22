set_path <- function(file, add.folder = NULL) {
        f_sep <- .Platform$file.sep
        file_name <- unlist(strsplit(file, f_sep))
        file_name <- file_name[length(file_name)]
        curr_wd <- unlist(strsplit(getwd(), f_sep))
        wdir <- grepl(" ", curr_wd)
        
        curr_wd[wdir] <-
                stringr::str_replace(
                        string      = curr_wd[wdir],
                        replacement = paste0("'", curr_wd[wdir], "'"),
                        pattern     = curr_wd[wdir]
                )
        
        curr_wd <- paste0(curr_wd, collapse = f_sep)
        
        if (!is.null(add.folder)) {
                add.folder <- paste0(f_sep, add.folder)
                return(paste0(curr_wd, add.folder, file_name))
                
        } else {
                return(paste0(curr_wd, file_name))
        }
}