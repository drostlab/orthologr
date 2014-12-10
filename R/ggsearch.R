#' @title Interface to the commandline tool ggsearch36.
#' @description This function is a wrapper for the tool ggsearch that performs
#' global alignment for each sequence in the inputfile with each sequence in the 
#' library and returns the alignment with the best score.
#' @param file a character string specifying the path to the sequence file.
#' @param library_file a character string specifying the path to the library file.
#' @param lib_format a character string specifying the file type of the library e.g. fasta.
#' @param options a character string specifying additional options for ggsearch
#' @param seq_type a character string specifying the type of the sequences
#' @param path a character string specifying the path to ggsearch if not included in 
#' the global environment
#' @param parse_output_to a charcater string specifying a file where the output of
#' ggsearch is parsed to
#' @details "Perform a rigorous Global/Global (GGSEARCH) [...] alignment between 
#' a protein sequence and another protein sequence or a protein database [or
#' a DNA sequence to a DNA sequence database ]. GGSEARCH implements a Needleman-Wunsch like alignment, except that affine gap 
#' penalties are used. GLSEARCH is most appropriate for global searches with a 
#' domain, which require local alignments within proteins. 
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @references ggsearch is part of the fasta36 package that can be found here 
#' http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml
#' @examples \dontrun{
#' # performing ggsearch on two protein files
#' ggsearch( file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'))
#'   
#' # performing ggsearch on two cds files parsing the output                   
#' ggsearch( file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           seqtype="cds", parse_output_to = "ggsearch.out")
#'               
#' }
#' @export
ggsearch <- function(file, library_file, options="", lib_format="fasta",
                     seq_type = "protein",
                     path=NULL, parse_output_to=NULL){
        
        if(seq_type == "protein"){
                options <- paste0(options, "-p")
        }
        
        operating_sys <- Sys.info()[1]
        
        if (operating_sys == "Darwin"){ 
                os <- "Mac"
                calc <- "ggsearch36"
        }
        
        if (operating_sys == "Linux"){
                os <- "Linux"
                calc <- "ggsearch36"
        }
        
        if (operating_sys == "Windows"){ 
                os <- "Windows"
                calc <- "ggsearch36.exe"
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