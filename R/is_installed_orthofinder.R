is_installed_orthofinder <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                sys_out <-
                        system("/opt/miniconda3/bin/orthofinder", intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have OrthoFinder2 installed locally on your machine or the PATH variable to the OrthoFinder2 program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "OrthoFinder")))
                return(TRUE)
        
}