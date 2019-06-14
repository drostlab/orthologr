is_installed_clustalo <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("clustalo --help", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have Clustal-Omega installed locally on your machine or the PATH variable to the Clustal-Omega program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "Clustal-Omega")))
                return(TRUE)
        
}