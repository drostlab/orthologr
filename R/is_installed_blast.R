is_installed_blast <- function(path = NULL) {
        # test if a valid BLAST version is installed
        tryCatch({
                if (is.null(path)) {
                        sys_out <-
                                system("blastp -version", intern = TRUE)
                } else {
                        sys_out <-
                                system(paste0(
                                        'export PATH=$PATH:',
                                        path, "'; blastp -version '"), intern = TRUE)
                }
                
                
        }, error = function(e)
                stop(
                        "It seems like you don't have BLAST installed locally on your machine or the PATH variable to the BLAST program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "blast")))
                return(TRUE)
        
}