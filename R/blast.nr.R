#' @title Perform a BLASTp search against NCBI nr
#' @description 
#' This function performs a BLASTp search of a given set of sequences NCBI nr.
#' @param query_file a character string specifying the path to the query file of interest (query organism).
#' @param nr.path path to the NCBI nr database folder.
#' @param seq_type a character string specifying the sequence type stored in the input file. 
#' Options are are: "cds", "protein", or "dna". 
#' In case of "cds", sequence are translated to protein sequences, in case of "dna", 
#' cds prediction is performed on the corresponding sequences which subsequently are 
#' translated to protein sequences. Default is \code{seq_type} = \code{"protein"}.
#' @param format a character string specifying the file format of the sequence file, 
#' e.g. \code{format} = \code{"fasta"}. Default is \code{format} = \code{"fasta"}.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is NULL, implicating that a set of default parameters is used when running BLAST.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to run BLAST searches.
#' @details
#' This function provides users an R interface to perform BLASTp searches against the NCBI nr database. 
#' @author Hajk-Georg Drost
#' @examples
#' 
#' \dontrun{
#' 
#' # perform a blastp serach of A. thaliana genes against NCBI nr
#'    blast.nr(
#'             query_file      = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'             nr.path         = "path/to/nr/database/folder",
#'             seq_type        = "protein",
#'             max.target.seqs = 1, 
#'             comp_cores      = 2 )
#'             
#' }
#' 
#' @references
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK1763/}
#' @seealso \code{\link{delta.blast}}, \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_blast}}, \code{\link{advanced_blast}}
#' @export


blast.nr <- function(query_file, 
                     nr.path        = NULL,
                     seq_type        = "protein", 
                     format          = "fasta",
                     eval            = "1E-5",
                     max.target.seqs = 500,
                     path            = NULL, 
                     blast_params    = NULL,
                     comp_cores      = 1){
        
        
        if(is.null(nr.path))
                stop("Please specify the path to your NCBI nr database.")
        
        nr.output <- data.table::copy(
                
                advanced_blast(
                        query_file      = query_file, 
                        subject_file    = "nr", 
                        seq_type        = seq_type, 
                        format          = format,
                        path            = path, 
                        blast_params    = paste0("-evalue ",eval," -num_threads ",
                                                 comp_cores," -max_target_seqs ",
                                                 max.target.seqs," ",blast_params),
                        blast_algorithm = "blastp",
                        db_path         = nr.path )
                
        )
        
        return(nr.output)
        
}








