is_installed_orthofinder <- function(path = NULL) {
        # test if a valid BLAST version is installed
        tryCatch({
                if (is.null(path)) {
                        sys_out <-
                                system("/opt/miniconda3/bin/orthofinder", intern = TRUE)
                } else {
                        sys_out <-
                                system(file.path(path, "orthofinder"), intern = TRUE)
                }
                
        }, error = function(e)
                stop(
                        "It seems like you don't have OrthoFinder2 installed locally on your machine or the PATH variable to the OrthoFinder2 program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "OrthoFinder")))
                return(TRUE)
        
}