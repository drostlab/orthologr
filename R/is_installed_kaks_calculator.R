is_installed_kaks_calculator <- function() {
        # test if a valid BLAST version is installed
        tryCatch({
                operating_sys <- Sys.info()[1]
                
                if (operating_sys == "Darwin") {
                        sys_out <-
                                system("KaKs_Calculator -h", intern = TRUE)
                }
                
                if (operating_sys == "Linux") {
                        sys_out <-
                                system("KaKs_Calculator -h", intern = TRUE)
                }
                
                if (operating_sys == "Windows") {
                        sys_out <-
                                system("KaKs_Calculator.exe -h", intern = TRUE)
                }
                
        }, error = function(e)
                stop(
                        "It seems like you don't have KaKs_Calculator installed locally on your machine or the PATH variable to the BLAST program is not set correctly. Please read the INSTALLATION vignette for details.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "Kaks_Calculator")))
                return(TRUE)
        
}