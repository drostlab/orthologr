#' @title Function to get a codon alignment.
#' @description This function takes a protein alignment and the corresponding coding sequences 
#' both given as files and uses the method specified by tool to create the corresponding codon
#' alignment.
#' @param file_aln a character string specifying the path to the file storing the protein alignment in CLUSTAL or FASTA format.
#' @param file_nuc a character string specifying the path to the file storing the coding sequences in multiple FASTA format.
#' @param format a character string specifying the file format used to store the codon alignment, e.g. "fasta", "clustal".
#' @param tool a character string specifying the program that should be used e.g. "pal2nal". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the codon alignment program (in case you don't use the default path).
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common codon alignment tools such as "PAL2NAL".
#' @return if get_aln is TRUE an object of class alignment of the seqinr package.
#' @export
codon_aln <- function(file_aln, file_nuc, format = "clustal", tool, get_aln = FALSE, path){
        
        # the Pal2Nal program is stored as executable within
        # the R package environment: 'exec' folder
        # this is not an elegant version, only motivated by this discussion:
        # http://stackoverflow.com/questions/13463103/inst-and-extdata-folders-in-r-packaging
        # is there a better choice?
        #path <- paste0(system.file(package = "orthologr"),"/exec/pal2nal.v14/")
        #path <- "/exec/pal2nal.v14/"
        
        if(!is.element(tool,c("pal2nal")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(format,c("clustal", "fasta")))
                stop("Please choose a format that is supported by this function.")
        
        if(!file.exists("_alignment/")){
                
                dir.create("_alignment")
        }
        
        # RIGHT NOW EACH NEW RUN OF THE FUNCTION OVERWRITES
        # THE EXISTING *.aln FILE
        # IN FUTURE VERSIONS WE SHOULD TRY TO KEEP (STORE) EXISTING FILES
        # AND RENAME THE NEW ONES
        
        if(tool == "pal2nal"){
        
                file.out <- "_alignment/pal2nal.aln"
                
                # test whether the connection to pal2nal works
                tryCatch(
                        
                        system(paste0(path,"pal2nal.pl ",file_aln," ",file_nuc," >",file.out))
                
                , stop( paste0("Please check the correct path to ",tool,
                               "... the interface call did not work properly.") )
                        
                )
                
                if(get_aln){
                        dna_aln <- seqinr::read.alignment(file = file.out, format = "clustal")
                        return(dna_aln)
                }
        }
}

