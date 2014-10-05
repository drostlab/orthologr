#' @title An interface function to common pairwise alignment tools.
#' @description This function takes a FASTA file containing two DNA or amino acid sequences
#' that shall be aligned and computes a paiwise alignment using a defined alignment method.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program/algorithm that should be used: "NW".
#' @param seq_type a character string specifying the sequence type stored within the given FASTA file. Options are
#' \code{seq_type} = "protein", "cds", "dna". Default is \code{seq_type} = "protein". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param paiwise_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{pairwise_aln_name} = \code{NULL} denoting a default name: 'toolname_seq_type.aln' .
#' @param path a character string specifying the path to the pairwise alignment program (in case you don't use the default path).
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details This function provides an interface between R and common pairwise alignment programs.
#' @examples \dontrun{        
#'                                             
#' # Needleman-Wunsch Example:  
#' 
#' # in case Biostrings works properly
#' pairwise_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'              tool = "NW", get_aln = TRUE, seq_type = "protein")
#' 
#'                                                    
#'                                                                                                                                                          
#' }
#' @return In case the argument \code{get_aln} is set \code{TRUE}, an object of class alignment of the seqinr package is returned.
#' @export
pairwise_aln <- function(file, tool = "NW", seq_type = "protein", 
                         get_aln = FALSE, pairwise_aln_name = NULL,
                         path = NULL, quiet = FALSE){
        
        
        if(!is.pairwise_aln_tool(tool))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(seq_type,c("dna","cds","protein")))
                stop("Please choose a seq_type that is supported by this function.")
        
        if(is.element(seq_type,c("dna","cds")))
                seqtype <- "DNA"
           
        if(seq_type == "protein")
                seqtype <- "AA"
        
        
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
                if(is.element(seq_type,c("dna","cds"))){
#                         not possible as is sorts the input by name
#                         seqs <- read.cds(file=file, format=format)
#                         names <- dt[,geneids] 
#                         seqs <- sapply(dt[,seqs] , Biostrings::DNAString)
                        
                        input <- seqinr::read.fasta(file=file, seqtype = seqtype)
                        names <- names(input) 
                        seqs <- lapply(input, function(x){return (Biostrings::DNAString(seqinr::c2s(x)))})      
                }
                if(seq_type == "protein"){
#                         dt <- read.proteome(file=file, format=format)
#                         names <- dt[,geneids] 
#                         seqs <- sapply(dt[,seqs] , Biostrings::AAString)

                        input <- seqinr::read.fasta(file=file, seqtype = seqtype)
                        names <- names(input) 
                        seqs <- lapply(input, function(x){return (Biostrings::AAString(seqinr::c2s(x)))})      
                }
                
                # align
                aln  <- Biostrings::pairwiseAlignment(pattern = seqs[[1]], subject = seqs[[2]],
                                                      type = "global")
                 
                #write file -> comment Hajk: what does pattern() and subject() do in pattern(aln), subject(aln) ?
                # pattern(aln) == get aligned sequence that was set as pattern input
                # subject(aln) == get aligned sequence that was set as subject input
                seqinr::write.fasta(sequences = list(pattern(aln), subject(aln)) , names = names,
                                    file.out = file.out)
                
                if(!quiet){ print(paste0("File successfully written to ",file.out)) }
                
                if(get_aln) return(aln)
        }
}




