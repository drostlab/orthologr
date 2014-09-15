#' @title Function to get a multiple alignment.
#' @description This function takes a multiple FASTA file containing DNA or aminoacid sequences
#' to be aligned and proceeds the method defined by tool to align them.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param file.out a character string specifying the path where the alignment should be stored.
#' @param tool a character string specifying the program that should be used e.g. "clustalw". 
#' @param get.aln a logical value indicating whether the produced alignment should be returned.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common multiple alignment programs
#' such as "clustalw", "tcoffee", "muscle", and "clustalo".
#' @return if get.aln is TRUE an object of class alignment of the seqinr package.
#' @export
multi_aln<-function(file, tool, get.aln="FALSE"){
        
        if(!is.element(tool,c("clustalw", "tcoffee", "muscle", "clustalo")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!file.exists("_alignment/")){
                
                dir.create("_alignment")
        }
        
        if(tool=="clustalw"){
        
                file.out <- "_alignment/clustalw.aln"
                
                system(paste0("clustalw ",file," -outfile=",file.out," -quiet"))

                if(get.aln){
                        aln <- seqinr::read.alignment(file.out, format="clustal")
                        return(aln)
                }
                
                # print(paste0("Alignment saved in ", file.out))
        }
 
        if(tool=="tcoffee"){
                
                file.out <- "_alignment/tcoffee.aln"
                
                system(paste0("t_coffee ",file," >",file.out))
                
                if(get.aln){
                        aln <- seqinr::read.alignment(file.out, format="clustal")
                        return(aln)
                }
                
                # print(paste0("Alignment saved in ", file.out))
        }
 
         if(tool=="muscle"){
                 
                 file.out <- "_alignment/muscle.aln"
                 
                 system(paste0("muscle -in ",file," -out ",file.out))
                 
                 if(get.aln){
                         aln <- seqinr::read.alignment(file.out, format="fasta")
                         return(aln)
                 }
                 
                 # print(paste0("Alignment saved in ", file.out))
         }
 
         if(tool=="clustalo"){
                
                 file.out <- "_alignment/clustalo.aln"
                 
                 system(paste0("clustalo -i ",file," -o ",file.out," --outfmt clustal"))
                 
                 if(get.aln){
                         aln <- seqinr::read.alignment(file.out, format="clustal")
                         return(aln)
                 }
                 
                 # print(paste0("Alignment saved in ", file.out))
         }

}




