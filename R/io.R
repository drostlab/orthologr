#' @title Read the genome of a given organism
#' @description This function reads an organism specific genome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the genome.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readDNAStringSet}} function.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details The \code{read.genome} function takes a string specifying the path to the genome file
#' of interest as first argument.
#'
#' It is possible to read in different genome file standards such as \emph{fasta} or \emph{fastq}.
#' Genomes stored in fasta files can be downloaded from http://ensemblgenomes.org/info/genomes.
#'
#' @examples \dontrun{
#' # reading a genome stored in a fasta file
#' Ath.genome <- read.genome(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                            format = "fasta")
#' }
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.genome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","fastq")))
                stop("Please choose a file format that is supported by this function.")
        
        geneids <- seqs <- NULL
        
        if(format == "fasta"){
                
                tryCatch({
                        
                        genome <- Biostrings::readDNAStringSet(filepath = file, format = format, ...)
                        genome_names <- as.vector(unlist(sapply(genome@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
                        genome.dt <- data.table::data.table(geneids = genome_names,
                                                            seqs = toupper(as.character(genome)))
                        
                        data.table::setkey(genome.dt,geneids)
                        
                }, error = function(e){ stop(paste0("File ",file, " could not be read properly. \n",
                                                    "Please make sure that ",file," contains only DNA sequences and is in ",format," format."))}
                )
        }
        
        return(genome.dt)
}


#' @title Read the proteome of a given organism
#' @description This function reads an organism specific proteome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the proteome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readAAStringSet}} function.
#' @author Hajk-Georg Drost
#' @details The \code{read.proteome} function takes a string specifying the path to the proteome file
#' of interest as first argument.
#'
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{fastq}.
#'
#' Proteomes stored in fasta files can be downloaded from http://www.ebi.ac.uk/reference_proteomes.
#'
#' @examples \dontrun{
#' # reading a proteome stored in a fasta file
#' Ath.proteome <- read.proteome(system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'                                format = "fasta")
#' }
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.proteome <- function(file, format, ...){
        
        if(!is.element(format,c("fasta","fastq")))
                stop("Please choose a file format that is supported by this function.")
        
        geneids <- seqs <- NULL 
        
                
                tryCatch({
                        
                        proteome <- Biostrings::readAAStringSet(filepath = file, format = format, ...)
                        proteome_names <- as.vector(unlist(sapply(proteome@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
                        proteome.dt <- data.table::data.table(geneids = proteome_names,
                                                              seqs = toupper(as.character(proteome)))
                        
                        data.table::setkey(proteome.dt, geneids)
                        
                }, error = function(e){ stop("File ",file, " could not be read properly.", "\n",
                                                    "Please make sure that ",file," contains only amino acid sequences and is in ",format," format.")}
                )
        
        return(proteome.dt)
}


#' @title Read the CDS of a given organism
#' @description This function reads an organism specific CDS stored in a defined file format.
#' @param file a character string specifying the path to the file storing the CDS.
#' @param format a character string specifying the file format used to store the CDS, e.g. "fasta", "fatsq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readDNAStringSet}} function.
#' @author Hajk-Georg Drost
#' @details The \code{read.cds} function takes a string specifying the path to the cds file
#' of interest as first argument.
#'
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{fastq}.
#'
#' CDS stored in fasta files can be downloaded from http://www.ensembl.org/info/data/ftp/index.html.
#'
#' @examples \dontrun{
#' # reading a cds file stored in fasta format
#' Ath.cds <- read.cds(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                     format = "fasta")
#' }
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

read.cds <- function(file, format, delete_corrupt = TRUE, ...){
        
        if(!is.element(format,c("fasta","fastq")))
                stop("Please choose a file format that is supported by this function.")
        
        
        geneids <- seqs <- NULL
        
                tryCatch({
                        
                        cds_file <- Biostrings::readDNAStringSet(filepath = file, format = format, ...)
                        
                        cds_names <- as.vector(unlist(sapply(cds_file@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
                        
                        cds.dt <- data.table::data.table(geneids = cds_names ,
                                                         seqs = tolower(as.character(cds_file)))
                        
                        
                        data.table::setkey(cds.dt, geneids)
                        
                        mod3 <- function(x) { return((nchar(x) %% 3) == 0) }
                        
                        all_triplets <- cds.dt[ , mod3(seqs)]
                        
                                           
                }, error = function(e){ stop("File ",file, " could not be read properly.","\n",
                                                    "Please make sure that ",file," contains only CDS sequences and is in ",format," format.")}
                )
        
        if(!all(all_triplets)){
                
                if(delete_corrupt)
                        return(cds.dt[-which(!all_triplets) , list(geneids, seqs)])
                
                if(!delete_corrupt)
                        return(cds.dt)
                
        } else {
                
                return(cds.dt)
        }
        
        
}

#' @title Translate DNA to Amino Acids
#' @description This function takes CDS sequence as string as input an returns the
#' corresponding amino acid sequence as string. 
#' @param sequence a character string specifying CDS sequence of interest.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @return A character string specifying the corresponding amino acid sequence.
#' @references \code{\link[seqinr]{translate}}
#' @examples
#' 
#' # an example DNA sequence
#' DNA <- c("ACCGGTTTAAAGGCGTTA")
#' 
#' # translating DNA to a protein sequence
#' transl(DNA)
#' 
#' @export
transl <- function(sequence){
        return(seqinr::c2s(seqinr::translate(seqinr::s2c(sequence))))
}

