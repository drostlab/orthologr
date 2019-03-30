is_installed_blast <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("blastp -version", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have BLAST installed locally on your machine or the PATH variable to the BLAST program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "blast")))
                return(TRUE)
        
}