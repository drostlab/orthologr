
#' @title Function to read a genome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details The \code{read.genome} function takes a string specifying the path to the genome file
#' of interest as first argument.
#' 
#' It is possible to read in different genome file standards such as \emph{fasta} or \emph{genebank}.
#' Genomes stored in fasta files can be downloaded from http://ensemblgenomes.org/info/genomes.
#' 
#' @examples \dontrun{
#' # reading a genome stored in a fasta file
#' Aly.genome <- read.genome("ortho_lyra_cds.fasta", format = "fasta")
#' }
#' 
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export 
read.genome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                genome <- vector(mode = "list")
                genome <- seqinr::read.fasta(file, seqtype = "DNA", ...)
                genome.dt <- data.table::data.table(geneids = names(genome), 
                                                      seqs = lapply(genome, seqinr::c2s))
                data.table::setkey(genome.dt,geneids)
        }

        return(genome.dt)
}



#' @title Function to read a proteome of a given organism
#' @description This function reads an organism specific proteome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the proteome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details The \code{read.proteome} function takes a string specifying the path to the proteome file
#' of interest as first argument.
#' 
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{genebank}.
#' 
#' Proteomes stored in fasta files can be downloaded from http://www.ebi.ac.uk/reference_proteomes.
#' 
#' @examples \dontrun{
#' # reading a proteome stored in a fasta file
#' Aly.proteome <- read.proteome("ortho_lyra_aa.fasta", format = "fasta")
#' }
#' 
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export 

read.proteome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                proteome <- vector(mode = "list")
                proteome <- seqinr::read.fasta(file, seqtype = "AA", ...)
                proteome.dt <- data.table::data.table(geneids = names(proteome), 
                                       seqs = unlist(lapply(proteome, seqinr::c2s)))
                data.table::setkey(proteome.dt,geneids)
        }
            
        return(proteome.dt)
}



#' @title Function to read the CDS of a given organism
#' @description This function reads an organism specific CDS stored in a defined file format.
#' @param file a character string specifying the path to the file storing the CDS.
#' @param format a character string specifying the file format used to store the CDS, e.g. "fasta", "gbk".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details The \code{read.cds} function takes a string specifying the path to the cds file
#' of interest as first argument.
#' 
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{genebank}.
#' 
#' CDS stored in fasta files can be downloaded from http://www.ensembl.org/info/data/ftp/index.html.
#' 
#' @examples \dontrun{
#' # reading a cds file stored in fasta format
#' Aly.cds <- read.cds("ortho_lyra_cds.fasta", format = "fasta")
#' }
#' 
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export 
read.cds <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","gbk")))
                stop("Please choose a file format that is supported by this function.")
        
        
        if(format == "fasta"){
                cds <- vector(mode = "list")
                cds <- seqinr::read.fasta(file, seqtype = "DNA", ...)
                cds.dt <- data.table::data.table(geneids = names(cds), 
                                                      seqs = unlist(lapply(cds, seqinr::c2s)))
                data.table::setkey(cds.dt,geneids)
        }
               
        return(cds.dt)
}


#' @title Interface function to BLAST+
#' @description This function performs a blast+ search of a given set of protein sequences against a given database.
#' @param queries 
#' @param database 
#' @param eval 
#' @param blast_name
#' @param blast_out
#' @param path
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @examples \dontrun{
#' # performing a best hit BLAST search
#' blast( set("ortho_lyra_cds.fasta") , set("ortho_thal_cds.fasta", makedb = TRUE))
#' }
#' 
#' @return A data.table storing the geneids of orthologous genes (best hit) 
#' in the first column and the amino acid sequences in the second column.
#' @export 
blast <- function(queries, database, eval = "1E-5",
                  blast_name = "_blast/_blastinput.fasta",
                  blast_out = "_blast/_blastresult.csv",
                  path = NULL){
        
        nrows <- nrow(queries)
        
        res <-NULL
        if(!file.exists("_blast/")){ 
                
                dir.create("_blast")
        
        }
        
        # here to come: a parallelized version of the BLAST process
        testAA <- as.list(queries[ , 3])
        name <- sapply(queries[ , 1], FUN = toString)
        seqinr::write.fasta(testAA, names = name,
                            nbchar = 80,open = "w", 
                             file.out = blast_name)
                        if(is.null(path)){
                                
                                system( 
                                        paste0("blastp -db ",database," -query ",blast_name,
                                               " -evalue ",eval," -out ", blast_out ," -outfmt 6")
                                )
                        
                        } else {
                                
                                system( 
                                        paste0("export PATH=$PATH:",path,"; blastp -db ",
                                               database," -query ",blast_name," -evalue ",
                                               eval," -out ", blast_out ," -outfmt 6")              
                                )
                        
                        }
        

                hit_table <- data.table::fread(input = blast_out, sep = "\t", header = FALSE)
                #setnames(hit_table, new = c())
                #setkey(hit_table, geneids)
                
                uniques <- vector(mode = "character") 
                uniques <- unique(hit_table[ , 1])
                
                for( entry in uniques ){
                        res <-rbind(res, hit_table[hit_table[,1]==entry,][1,1:2] )
                }

        return(res)
}

#' @title Function to perform a BLAST best hit search
#' @description This function performs a blast+ search (best hit) of a given set of protein sequences against a given database.
#' @param A an object returned by the set() function.
#' @param B an object returned by the set() function having the makedb parameter set to TRUE.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details Given a set of protein sequences A, a best hit blast search is being performed from A to database.
#' @examples \dontrun{
#' # performing gene orthology inference using the best hit (BH) method 
#' blast_best( set("ortho_lyra_aa.fasta") , set("ortho_thal_aa.fasta", makedb = TRUE))
#' }
#' 
#' @return A data.table as returned by the blast() function, storing the geneids 
#' of orthologous genes (best hit) in the first column and the amino acid sequences in the second column.
#' @export 
blast_best <- function(A,B){
        
        return(blast(queries = A[[1]], database = B[[2]]))
        
}


#' @title Function to perform a BLAST reciprocal best hit (RBH) search
#' @description This function performs a blast+ search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
#' @param A an object returned by the set() function.
#' @param B an object returned by the set() function having the makedb parameter set to TRUE.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details Given a set of protein sequences A and a different set of protein sequences B,
#'  first a best hit blast search is being performed from A to B: blast(A,B) and afterwards
#'  a best hit blast search is being performed from B to A: blast(B,A). Only protein sequences 
#'  that were found to be best hits in both directions are retained and returned.
#'  
#'  
#' @examples \dontrun{
#' # performing gene orthology inference using the best hit (BH) method 
#' blast_rec( set("ortho_lyra_aa.fasta") , set("ortho_thal_aa.fasta", makedb = TRUE))
#' }
#' 
#' @return A data.table as returned by the blast() function, storing the geneids 
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @export 
blast_rec <- function(A, B){
        orthoA <- blast_best(A,B)
        orthoB <- blast_best(B,A)
        return ( merge(orthoA, orthoB, by.x = c(1,2), by.y = c(2,1)))
}



