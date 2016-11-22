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