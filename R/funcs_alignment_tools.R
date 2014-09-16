#' @title Function to get a multiple alignment.
#' @description This function takes a multiple FASTA file containing DNA or aminoacid sequences
#' to be aligned and proceeds the method defined by tool to align them.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program that should be used e.g. "clustalw". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common multiple alignment programs
#' such as "clustalw", "tcoffee", "muscle", and "clustalo".
#' @return If get_aln is TRUE an object of class alignment of the seqinr package.
#' @export
multi_aln <- function(file, tool, get_aln = "FALSE", path = NULL){
        
        if(!is.element(tool,c("clustalw", "tcoffee", "muscle", "clustalo")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!file.exists("_alignment/")){
                
                dir.create("_alignment")
        }
        
        # RIGHT NOW EACH NEW RUN OF THE FUNCTION OVERWRITES
        # THE EXISTING *.aln FILE
        # IN FUTURE VERSIONS WE SHOULD TRY TO KEEP (STORE) EXISTING FILES
        # AND RENAME THE NEW ONES
        
        if(tool == "clustalw"){
        
                file.out <- "_alignment/clustalw.aln"
                
                # test whether the connection to clustalw works
                tryCatch(
                
                       if(is.null(path)){
                        
                                system(paste0("clustalw ",file," -outfile=",file.out," -quiet"))
                        
                        } else {
                                system(
                                        paste0("export PATH=$PATH:",path,"; ","clustalw ",
                                               file," -outfile=",file.out," -quiet")
                                )
                        }
                , stop( paste0("Please check the correct path to ",tool,
                              "... the interface call did not work properly.") )
                       
                       )
                
                if(get_aln){
                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                        return(aln)
                }
        }
 
        if(tool == "tcoffee"){
                
                file.out <- "_alignment/tcoffee.aln"
                
                # test whether the connection to t_coffee works
                tryCatch(
                
                        if(is.null(path)){
                        
                                system(paste0("t_coffee ",file," >",file.out))
                        
                        } else {
                        
                                system(paste0("export PATH=$PATH:",path,"; ","t_coffee ",
                                              file," >",file.out))
                        
                        }
                
                , stop( paste0("Please check the correct path to ",tool,
                               "... the interface call did not work properly.") )
                
                )
                
                if(get_aln){
                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                        return(aln)
                }
        }
 
         if(tool == "muscle"){
                 
                 file.out <- "_alignment/muscle.aln"
                 
                 # test whether the connection to muscle works
                 tryCatch(
                         
                          if(is.null(path)){
                         
                                   system(paste0("muscle -in ",file," -out ",file.out))
                 
                           } else {
                         
                                   system(paste0("export PATH=$PATH:",path,"; ","muscle -in ",
                                       file," -out ",file.out))
                         
                           }
                 
                 , stop( paste0("Please check the correct path to ",tool,
                                "... the interface call did not work properly.") )
                 
                 )
                 
                 if(get_aln){
                         aln <- seqinr::read.alignment(file.out, format = "fasta")
                         return(aln)
                 }
         }
 
         if(tool == "clustalo"){
                
                 file.out <- "_alignment/clustalo.aln"
                 
                 
                 # test whether the connection to clustalo works
                 tryCatch(
                         
                         if(is.null(path)){
                         
                                 system(paste0("clustalo -i ",file," -o ",file.out," --outfmt clustal"))
                 
                         } else {
                         
                                 system(paste0("export PATH=$PATH:",path,"; ","clustalo -i ",
                                               file," -o ",file.out," --outfmt clustal"))
                         
                        }
                 
                 , stop( paste0("Please check the correct path to ",tool,
                                "... the interface call did not work properly.") )
                 
                 )
                 
                 if(get_aln){
                         aln <- seqinr::read.alignment(file.out, format = "clustal")
                         return(aln)
                 }
         }

}




