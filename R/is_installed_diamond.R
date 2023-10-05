# from rdiamond (github.com/drostlab/rdiamond)
is_installed_diamond <- function(diamond_exec_path = NULL) {
        # test if a valid DIAMOND version is installed
        tryCatch({
                if (!is.null(diamond_exec_path))
                        # e.g. diamond_exec_path = '/opt/miniconda3/bin/'
                        sys_out <-
                                system(paste0(file.path(diamond_exec_path,"diamond version")), intern = TRUE)
                
                if (is.null(diamond_exec_path))
                        # e.g. diamond_exec_path = '/opt/miniconda3/bin/'
                        sys_out <-
                                system("diamond version", intern = TRUE)
                
        }, error = function(e)
                stop(
                        "It seems like you don't have DIAMOND2 installed locally on your machine or the PATH variable to the DIAMOND2 program is not set correctly.",
                        " Please use the 'diamond_exec_path' argument to specify a path to your diamond excecutable, e.g. diamond_exec_path = '/usr/local/bin' or diamond_exec_path = '/opt/miniconda3/bin'.",
                        call. = FALSE
                ))
        
        if (any(stringr::str_detect(sys_out, "diamond")))
                return(TRUE)
        
}