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
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
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

read.cds <- function(file, format, delete_corrupt_cds = TRUE, ...) {
        if (!is.element(format, c("fasta", "fastq")))
                stop("Please choose a file format that is supported by this function.")
        
        
        geneids <- seqs <- NULL
        
        tryCatch({
                cds_file <-
                        Biostrings::readDNAStringSet(filepath = file, format = format, ...)
                
                cds_names <-
                        as.vector(unlist(sapply(cds_file@ranges@NAMES, function(x) {
                                return(strsplit(x, " ")[[1]][1])
                        })))
                
                cds.dt <-
                        data.table::data.table(geneids = cds_names ,
                                               seqs = tolower(as.character(cds_file)))
                
                
                data.table::setkey(cds.dt, geneids)
                
                mod3 <-
                        function(x) {
                                return((nchar(x) %% 3) == 0)
                        }
                
                all_triplets <- cds.dt[, mod3(seqs)]
                
                n_seqs <- nrow(cds.dt)
                
        }, error = function(e) {
                stop(
                        "File ",
                        file,
                        " could not be read properly.",
                        "\n",
                        "Please make sure that ",
                        file,
                        " contains only CDS sequences and is in ",
                        format,
                        " format."
                )
        })
        
        if (!all(all_triplets)) {
                message(
                        "There seem to be ",
                        length(which(!all_triplets)),
                        " coding sequences in your input dataset which cannot be properly divided in base triplets, because their sequence length cannot be divided by 3."
                )
                corrupted_file <-
                        paste0(basename(file), "_corrupted_cds_seqs.fasta")
                
                message(
                        "A fasta file storing all corrupted coding sequences for inspection was generated and stored at '",
                        file.path(getwd(), corrupted_file),
                        "'."
                )
                message("\n")
                corrupted_seqs <- as.data.frame(cds.dt[which(!all_triplets)])
                seq_vector <- corrupted_seqs$seqs
                names(seq_vector) <- corrupted_seqs$geneids
                corrupted_seqs_biostrings <- Biostrings::DNAStringSet(seq_vector, use.names = TRUE)                                       
                Biostrings::writeXStringSet(corrupted_seqs_biostrings, filepath = corrupted_file)
                
                if (delete_corrupt_cds) {
                        message(
                                "You chose option 'delete_corrupt_cds = TRUE', thus corrupted coding sequences were removed.",
                                "If after consulting the file '",
                                corrupted_file,
                                "' you still wish to retain all coding sequences please specify the argument 'delete_corrupt_cds = FALSE'."
                        )
                        message("\n")
                        return(cds.dt[-which(!all_triplets) , list(geneids, seqs)])
                }
                
                if (!delete_corrupt_cds) {
                        message(
                                "You chose option 'delete_corrupt_cds = FALSE', thus corrupted coding sequences were retained for subsequent analyses.")
                        message(
                                "The following modifications were made to the CDS sequences that were not divisible by 3:")
                        message(
                                "- If the sequence had 1 residue nucleotide then the last nucleotide of the sequence was removed.")
                        message(
                                "- If the sequence had 2 residue nucleotides then the last two nucleotides of the sequence were removed.")
                        message(
                                "If after consulting the file '",
                                corrupted_file,
                                "' you wish to remove all corrupted coding sequences please specify the argument 'delete_corrupt_cds = TRUE'."
                        )
                        
                        mod3_residue_1 <-
                                function(x) {
                                        return((nchar(x) %% 3) == 1)
                                }
                        
                        mod3_residue_2 <-
                                function(x) {
                                        return((nchar(x) %% 3) == 2)
                                }
                        
                        residue_1 <- cds.dt[ , mod3_residue_1(seqs)]
                        residue_2 <- cds.dt[ , mod3_residue_2(seqs)]
                        
                        residue_1_seqs <- as.character(cds.dt[which(residue_1) , seqs])
                        residue_2_seqs <- as.character(cds.dt[which(residue_2) , seqs])
                        
                        residue_1_seqs_vec <- as.character(sapply(residue_1_seqs, function(x) {
                                stringr::str_sub(x, 1, nchar(x) - 1)
                        }))
                
                        residue_2_seqs_vec <- as.character(sapply(residue_2_seqs, function(x) {
                                stringr::str_sub(x, 1, nchar(x) - 2)
                        }))
                        
                        cds.dt[which(residue_1) , seqs := residue_1_seqs_vec]
                        cds.dt[which(residue_2) , seqs := residue_2_seqs_vec]
                        
                        all_triplets_new <- cds.dt[ , mod3(seqs)]
                        
                        if (any(!all_triplets_new)) {
                                stop("Something went wring during the trimming process. Not all sequences were trimmed properly.", call. = FALSE)
                        } else {
                                message("All corrupted CDS were trimmed.")
                        }
                        
                        n_seqs_new <- nrow(cds.dt)
                        
                        if(!(n_seqs == n_seqs_new))
                                stop("After trimming corrupted CDS some sequences seem to be lost. Please check what might have gone wrong with the sequence trimming.", call. = FALSE)
                        
                        return(cds.dt)
                }
                
        } else {
                return(cds.dt)
        }
        
}





