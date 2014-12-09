#' @title Interface to the commandline tool ssearch36.
#' @description This function is a wrapper for the tool ssearch that performs
#' global alignment for each sequence in the inputfile with each sequence in the 
#' library and returns the alignment with the best score.
#' @param file a character string specifying the path to the sequence file.
#' @param library_file a character string specifying the path to the library file.
#' @param lib_format a character string specifying the file type of the library e.g. fasta.
#' @param options a character string specifying additional options for ssearch
#' @param seq_type a character string specifying the type of the sequences
#' @param path a character string specifying the path to ssearch if not included in 
#' the global environment
#' @param parse_output_to a charcater string specifying a file where the output of
#' ssearch is parsed to
#' @details "Perform a rigorous Global/Global (ssearch) [...] alignment between 
#' a protein sequence and another protein sequence or a protein database [or
#' a DNA sequence to a DNA sequence database ]". ssearch is "using the Smith-Waterman
#' algorithm. ssearch36 uses SSE2 acceleration, and is only 2 - 5X slower than
#' fasta36
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @references ssearch is part of the fasta36 package that can be found here 
#' http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml
#' @examples 
#' # performing ssearch on two protein files
#' ssearch( file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'))
#'   
#' # performing ssearch on two cds files parsing the output                   
#' ssearch( file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           seqtype="cds", parse_output_to = "ssearch.out")
#'               
ssearch <- function(file, library_file, options="", lib_format="fasta",
                     seq_type = "protein",
                     path=NULL, parse_output_to=NULL){
        
        if(seq_type == "protein"){
                options <- paste0(options, "-p")
        }
        
        operating_sys <- Sys.info()[1]
        
        if (operating_sys == "Darwin"){ 
                os <- "Mac"
                calc <- "ssearch36"
        }
        
        if (operating_sys == "Linux"){
                os <- "Linux"
                calc <- "ssearch36"
        }
        
        if (operating_sys == "Windows"){ 
                os <- "Windows"
                calc <- "ssearch36.exe"
        }
        
        
        # differentiate system call
        
        if(is.null(path)){
                if(is.null(parse_output_to)){
                        system(paste0(calc," ",options," ",file," ",library_file))
                }
                if(!is.null(parse_output_to)){
                        system(paste0(calc," ",options," ",file," ",library_file," >",parse_output_to))
                }
        }
        if(!is.null(path)){
                if(is.null(parse_output_to)){
                        system(paste0(path,calc," ",options," ",file," ",library_file))
                }        
                if(!is.null(parse_output_to)){
                        system(paste0(path,calc," ",options," ",file," ",library_file," >",parse_output_to))
                }
        }
}