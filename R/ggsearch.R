#' @title Compute Global Alignments with GGSearch 
#' @description This function allows you to perform global alignments between pairs
#' of protein sequences or alignments between a protein sequence against a 
#' protein database based on the GGSeach program implemented in FASTA framework. 
#' 
#' @param file a character string specifying the path to the query sequence file.
#' @param library_file a character string specifying the path to the subject library file.
#' @param lib_format a character string specifying the file type of the library e.g. \code{lib_format} = \code{"fasta"}.
#' @param options a character string specifying additional options for GGSearch.
#' @param seq_type a character string specifying the type of biological sequences.
#' @param path a character string specifying the path to GGSearch. This argument only needs
#' to be specified when GGSearch cannot be executed from the default execution \code{PATH} (global environment).
#' @param parse_output_to a charcater string specifying a filename in which the output of
#' GGSearch is written.
#' @details 
#' 
#' This function takes two fasta files representing query and subject (library) sequences and computes
#' Global Pairwise Alignments.
#' 
#' As a result, this function returns the alignment score and additional parameters
#' computed by GGSearch. GGSearch implements a Needleman-Wunsch like alignment, except that affine gap 
#' penalties are used. 
#' 
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @references 
#' 
#' GGSearch is part of the fasta36 package that can be downloaded here 
#' 
#' \url{http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml}
#' 
#' Pearson WR and Lipman DJ. 1988. Improved tools for biological sequence comparison. PNAS 85: 2444-2448.
#' 
#'  \url{http://www.pnas.org/content/85/8/2444.long}
#' @examples \dontrun{
#' 
#' # performing GGSearch on two protein files
#' ggsearch( file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'))
#'   
#' # performing GGSearch on two cds files parsing the output                   
#' ggsearch( file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           library_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           seq_type = "cds", parse_output_to = "ggsearch.out")
#'               
#' }
#' @export

ggsearch <- function(file, library_file, options = "", lib_format = "fasta",
                     seq_type = "protein", path = NULL, parse_output_to = NULL){
        
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






