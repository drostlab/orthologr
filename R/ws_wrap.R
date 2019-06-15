#' @title Wrap whitespace in paths
#' @description Helper function for wrapping white spaces in folder paths.
#' @param paths a string or vector of strings containing paths.
#' @author Hajk-Georg Drost
#' @noRd

ws_wrap <- function(paths){
        
        paths <- stringr::str_split(paths,.Platform$file.sep)
        
        ws.paths <- unlist(lapply(paths, function(x) {
                
                stringr::str_c(unlist(sapply(x, function(word) {
                        if (stringr::str_detect(word," ")) 
                                return(stringr::str_replace(word,word,paste0("'",word,"'"))) else 
                                        return(word)})),collapse = .Platform$file.sep)
                
        }))
        
        return(ws.paths) 
}
