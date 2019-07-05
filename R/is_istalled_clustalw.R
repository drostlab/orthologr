is_installed_clustalw <- function() {
        # test if a valid BLAST version is installed
        operating_sys <- Sys.info()[1]
        
        if (operating_sys != "Windows") {
                call_clustalw <- "clustalw2"
                
        }
        
        if (operating_sys == "Windows")
                call_clustalw <- "clustalw2.exe"
        
        tryCatch({
                sys_out <-
                        system(paste0(call_clustalw, " -help"), intern = TRUE)
        }, error = function(e)
                stop(
                        "It seems like you don't have ClustalW installed locally on your machine or the PATH variable to the ClustalW program is not set correctly.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "CLUSTAL")))
                return(TRUE)
        
}