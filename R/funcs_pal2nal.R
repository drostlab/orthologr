#' @title Function to get a codon alignment.
#' @description This function takes a protein alignment and the corresponding coding sequences 
#' both given as files and uses the method specified by tool to create the corresponding codon
#' alignment.
#' @param file.aln a character string specifying the path to the file storing the protein alignment in CLUSTAL or FASTA format.
#' @param file.nuc a character string specifying the path to the file storing the coding sequences in multiple FASTA format.
#' @param file.out a character string specifying the path where the alignment should be stored.
#' @param format a character string specifying the file format used to store the codon alignment, e.g. "fasta", "clustal".
#' @param tool a character string specifying the program that should be used e.g. "pal2nal". 
#' @param get.aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifiying the absolute path where pal2nal.pl can be found.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common codon alignment tools such as "PAL2NAL".
#' @return if get.aln is TRUE an object of class alignment of the seqinr package.
#' @export
codon_aln<-function(file.aln, file.nuc, format="clustal", tool, get.aln=FALSE, path){
                
        if(!is.element(tool,c("pal2nal")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(format,c("clustal", "fasta")))
                stop("Please choose a format that is supported by this function.")
        
        if(!file.exists("_alignment/")){
                
                dir.create("_alignment")
        }
        
        
        if(tool=="pal2nal"){
        
                file.out <- "_alignment/pal2nal.aln"
                
                system(paste0(path,"pal2nal.pl ",file.aln," ",file.nuc," >",file.out))
                
                if(get.aln){
                        dna_aln <- seqinr::read.alignment(file=file.out, format="clustal")
                        return(dna_aln)
                }
        }
}

