is_installed_clustalw <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("clustalw -help", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have ClustalW installed locally on your machine or the PATH variable to the ClustalW program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "CLUSTAL")))
                return(TRUE)
        
}