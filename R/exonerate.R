#' @title Sequence Comparison with Exonerate
#' @description This function provides an interface to the Exonerate program: a generic tool for sequence alignment.
#' @param query_files a character vector specifying the paths to the fasta files of interest (query organism(s)).
#' @param subject_files a character vector specifying the paths to the fasta files of interest (subject organism(s)).
#' @param params a character string listing the input paramters that shall be passed to the executing Exonerate program. 
#' Default is \code{params} = \code{NULL}, implicating that a set of default parameters is used when running Exonerate.
#' @param core_interface a logical value enabling full flexibility to call exonerate with all additional parameters.
#' When \code{core_interface} = \code{TRUE}, only \code{system("exonerate ",params)} is being called and you can specify any
#' parameter constillation using the \code{params} argument. Note that you also have to specify the paths to you query and subject fasta files
#' when \code{core_interface} = \code{TRUE}. See examples for details.
#' @param format a character string specifying the file format of the sequence file. Right now only the 'fasta' format is supported.
#' @param clean_folders a boolean value specifying whether the internall "_exonerate" folder storing the output of exonerate shall be removed. 
#' Default is \code{clean_folders} = \code{TRUE}.
#' @details
#' 
#' When using this interface function it is assumed that you have the Exonerate program installed on your machine
#' and it can be executed from the default \emph{PATH} environment.
#' 
#' See \url{https://www.ebi.ac.uk/~guy/exonerate/exonerate.man.html} for detailed parameter specifications when using the Exonerate program.
#' 
#' 
#' Internally this function creates a folder "_exonerate" to store the corresponding output of
#' the Exonerate programs and will delete it (in case \code{clean_folders} = \code{TRUE}) after the output is
#' stored as R seqinr alignment object.
#' 
#' In case you want to retain the output of Exonerate, please select \code{clean_folders} = \code{FALSE}.
#' 
#' @author Hajk-Georg Drost
#' @return a sequence alignment as returned by the Exonerate program.
#' @references
#' 
#' Guy St C Slater and Ewan Birney (2005). Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6:31.
#' 
#' \url{https://www.ebi.ac.uk/~guy/exonerate/index.html}
#' 
#' @examples \dontrun{
#' 
#' # a simple ungapped alignment returned by exonerate
#' exonerate(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'))
#'           
#'           
#' ### peforming a Exonerate alignment using additional parameters
#' # 6-frame translated alignment: param = "--model coding2coding"
#' exonerate(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           params = "--model coding2coding")
#'
#' 
#' ### peforming a Exonerate with full flexibility
#' exonerate(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           params = paste0("--model coding2coding", 
#'           system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           " ",system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           core_interface = TRUE) 
#' 
#' }
#' @export

exonerate <- function(query_files, subject_files, params = NULL,core_interface = FALSE, format = "fasta", clean_folders = TRUE){
        
        if(!is.element(format,c("fasta")))
                stop("This function only supports the following input file formats: fasta")
        
        if(!file.exists("_exonerate"))
                dir.create("_exonerate")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        sessionName <- unlist(strsplit(query_files[1], f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
        sessionName <- sessionName[length(sessionName)]
        
        ExonerateOutputName <- paste0(sessionName,"_exonerate.fasta")
        
        tryCatch({
        
               if(is.null(params)){
                
                                    system(paste0("exonerate --query ",paste0(query_files,collapse = " "),
                                           " --target ", paste0(subject_files,collapse = " ")))
                }
        
                if(!is.null(params)){
                
                                    system(paste0("exonerate ",params," --query ",paste0(query_files,collapse = " "),
                                           " --target ", paste0(subject_files,collapse = " ")))
                }
        
        
                if(!is.null(params) & core_interface)
                                   system("exonerate ",params)
               
               
        
        },error = function(e){ clean_all_folders("_exonerate"); 
                               stop("The interface to the Exonerate program did not work properly. Is it executable from the default execution PATH?") }
        
        )
        
        # read exonerate output
        #aln <- seqinr::read.alignment(file = paste0("_exonerate",f_sep,ExonerateOutputName), format = format)
        
        if(clean_folders)
                clean_all_folders("_exonerate")
        
        
        return(aln)
        
}








