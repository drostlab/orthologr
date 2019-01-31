#' @title Helper function to select best BLAST hit based on minimum evalue
#' @details Extracts the best hit based on lowest e-value and longest splice variant 
#' when e-values are equal
#' @param x a tibble storing gene locus ids, splice variant ids, and blast output for filtering
#' @author Hajk-Georg Drost
filter_best_hits <- function(x) {
        
        min_val <- min(x$evalue)
        evalue <- alig_length <- NULL
        res <- dplyr::filter(x, evalue == min_val)
        if (nrow(res) > 1) {
                max_len <- max(res$alig_length)
                res <- dplyr::filter(res, alig_length == max_len)
        }
        
        if (nrow(res) > 1) 
                res <- dplyr::slice(res, 1)
        
        return(res)
}