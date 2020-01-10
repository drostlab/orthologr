is_installed_tcoffee <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("t_coffee -version", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have T-COFFEE installed locally on your machine or the PATH variable to the T-COFFEE program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "T-COFFEE")))
                return(TRUE)
        
}