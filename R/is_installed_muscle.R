is_installed_muscle <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("muscle -version", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have MUSCLE installed locally on your machine or the PATH variable to the MUSCLE program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "MUSCLE")))
                return(TRUE)
        
}