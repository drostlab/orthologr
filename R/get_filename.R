get_filename <- function(file_path){
        
        f_sep <- .Platform$file.sep
        split_path <- unlist(strsplit(file_path,f_sep))
        file_name <- unlist(strsplit(split_path[length(split_path)],"[.]"))[1]
        return(file_name)
}
