#' @title Save a proteome in fasta format
#' @description This function takes a proteome (in data.table notation) as returned by 
#' \code{\link{read.proteome}} and saves it to your hard drive.
#' @param proteome a data.table object storing amino acid sequences in the first column
#' and gene ids in the second column. See \code{\link{read.proteome}} for details.
#' @param file.name a character string specifying the name of the fasta output file.
#' @param nbchar number of characters per line.
#' @param open mode to open the output file, you can choose from "w" to write into a new file, or "a" to append at the end of an already existing file.
#' @param as.string when set to TRUE sequences are in the form of strings instead of vectors of single characters.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' \dontrun{
#' 
#' # read an example proteome
#' Ath.proteome <- read.proteome(
#'                       system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'), 
#'                       format = "fasta")
#' 
#' # use the write.proteome function to store it on your hard drive
#' write.proteome(Ath.proteome, "test_proteome.fasta")
#' 
#' }
#' 
#' @export

write.proteome <- function(proteome, file.name, nbchar = 80, open = "w", as.string = TRUE){
        # define visible bindings for global variables
        seqs <- geneids <- NULL
        seqinr::write.fasta( sequences = as.list(proteome[ , seqs]),
                             names     = proteome[ , geneids],
                             nbchar    = nbchar, 
                             open      = open,
                             file.out  = file.name,
                             as.string = as.string)
        
        
}