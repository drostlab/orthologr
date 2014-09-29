#' @title An interface function to common pairwise alignment tools.
#' @description This function takes a FASTA file containing two DNA or amino acid sequences
#' that shall be aligned and computes a paiwise alignment using a defined alignment method.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program/algorithm that should be used: "NW". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the pairwise alignment program (in case you don't use the default path).
#' @param paiwise_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{pairwise_aln_name} = \code{NULL} denoting a default name: 'toolname_seqtype.aln' .
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details This function provides an interface between R and common pairwise alignment programs.
#' @examples \dontrun{        
#'                                             
#' # Needleman-Wunsch Example:  
#' 
#' # in case Biostrings works properly
#' pairwise_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'              tool = "NW", get_aln = TRUE, seqtype="AA")
#'                                                    
#' }
#' @references 
#' @return In case the argument \code{get_aln} is set \code{TRUE}, an object of class alignment of the seqinr package is returned.
#' @export
pairwise_aln <- function(file, tool = "NW", seqtype, 
                         get_aln = FALSE, pairwise_aln_name = NULL,
                         path = NULL, quiet = FALSE){
        
        
        if(!is.pairwise_aln_tool(tool))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(seqtype,c("DNA","AA")))
                stop("Please choose a seqtype that is supported by this function.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!file.exists(paste0("_alignment",f_sep))){
                
                dir.create("_alignment")
        }
        
        if(!file.exists(paste0("_alignment",f_sep,"pairwise_aln",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"pairwise_aln"))
        }
        
        if(is.null(pairwise_aln_name)){
                
                file.out <- paste0("_alignment",f_sep,"pairwise_aln",f_sep,tool,"_",seqtype,".aln")
        }
        
        if(!is.null(pairwise_aln_name)){
                
                file.out <- paste0("_alignment",f_sep,"pairwise_aln",f_sep,pairwise_aln_name,"_",tool,"_",seqtype,".aln")    
        }
        
        # Needleman-Wunsch using Biostrings::pairwiseAlignment( type=global)
        if(tool == "NW"){
                
                #read file
                if(seqtype == "DNA"){
#                         not possible as is sorts the input by name
#                         seqs <- read.cds(file=file, format=format)
#                         names <- dt[,geneids] 
#                         seqs <- sapply(dt[,seqs] , Biostrings::DNAString)
                        
                        input <- seqinr::read.fasta(file=file, seqtype = seqtype)
                        names <- names(input) 
                        seqs <- lapply(input, function(x){return (Biostrings::DNAString(seqinr::c2s(x)))})      
                }
                if(seqtype == "AA"){
#                         dt <- read.proteome(file=file, format=format)
#                         names <- dt[,geneids] 
#                         seqs <- sapply(dt[,seqs] , Biostrings::AAString)

                        input <- seqinr::read.fasta(file=file, seqtype = seqtype)
                        names <- names(input) 
                        seqs <- lapply(input, function(x){return (Biostrings::AAString(seqinr::c2s(x)))})      
                }
                
                # align
                aln  <- Biostrings::pairwiseAlignment(pattern = seqs[[1]], subject = seqs[[2]],
                                                      type="global")
                 
                #write file
                seqinr::write.fasta(sequences = list(pattern(aln), subject(aln)) , names = names,
                                    file.out = file.out)
                
                if(!quiet){ print(paste0("File successfully written to ",file.out)) }
                
                if(get_aln) return(aln)
        }
}