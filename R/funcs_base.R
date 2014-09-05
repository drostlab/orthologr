
#' @title Function to check whether a given vector, data.frame, data.table, or list only stores DNA sequences
#' @description This function takes a vector, data.frame, data.table, or list
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @export 
check_DNA <- function(DNA_container){
        
        is.vector(DNA_container){
                
                
        }
        
        is.list(DNA_container){
                
                
                
        }
        
        
        is.data.frame(DNA_container){
                
                
                
        }
        
        
        is.data.table(DNA_container){
                
                
                
                
        }
        
}


#' @title Function to read a genome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @export 
read.genome <- function(file, format = "fasta", ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                genome <- vector(mode = "list")
                genome <- read.fasta(file, seqtype = "DNA", ...)
        }
        
        
        return(genome)
}



#' @title Function to read a proteome of a given organism
#' @description This function reads an organism specific proteome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the proteome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @export 

read.proteome <- function(file, format = "fasta", ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                proteome <- vector(mode = "list")
                proteome <- read.fasta(file, seqtype = "AA", ...)
        }
        
        
        return(proteome)
}



#' @title Function to read the CDS of a given organism
#' @description This function reads an organism specific CDS stored in a defined file format.
#' @param file a character string specifying the path to the file storing the CDS.
#' @param format a character string specifying the file format used to store the CDS, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @export 
read.cds <- function(file, format = "fasta", ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                cds <- vector(mode = "list")
                cds <- read.fasta(file, seqtype = "DNA", ...)
        }
        
        
        return(cds)
}





