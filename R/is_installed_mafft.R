is_installed_mafft <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("mafft -v", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have MAFFT installed locally on your machine or the PATH variable to the MAFFT program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "MAFFT")))
                return(TRUE)
        
}