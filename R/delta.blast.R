#' @title Perform a DELTA-BLAST Search
#' @description 
#' This function performs a DELTA-BLAST search of a given set of sequences against a given database
#' or against the cdd database itself.
#' @param query_file a character string specifying the path to the query file of interest (query organism).
#' @param subject_file a character string specifying the path to the subject file of interest (subject organism).
#' @param cdd.path path to the cdd database folder.
#' @param seq_type a character string specifying the sequence type stored in the input file. 
#' Options are are: "cds", "protein", or "dna". 
#' In case of "cds", sequence are translated to protein sequences, in case of "dna", 
#' cds prediction is performed on the corresponding sequences which subsequently are 
#' translated to protein sequences. Default is \code{seq_type} = \code{"cds"}.
#' @param format a character string specifying the file format of the sequence file, 
#' e.g. \code{format} = \code{"fasta"}. Default is \code{format} = \code{"fasta"}.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is NULL, implicating that a set of default parameters is used when running BLAST.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to run BLAST searches.
#' @param makedb_type a character string specifying the sequence type stored in the BLAST database that
#'  is generated using 'makeblastdb'. Options are: "protein" and "nucleotide". Default is \code{makedb_type} = \code{"protein"}.
#' @details
#' This function provides an interface between \code{delta-blast} and R. Given a query file, a delta-blast
#' search is perfored either against a subject database or against the cdd database itself.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # perform a delta-blast serach between A. thaliana and A. lyrata genes
#' delta.blast(
#'             query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'             cdd.path     = "path/to/cdd/database/folder" )
#'             
#' # perform a delta-blast serach between A. thaliana and the cdd database
#' delta.blast(
#'             query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'             subject_file = "cdd_delta",
#'             cdd.path     = "path/to/cdd/database/folder" )
#' 
#' @references
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.Performing_a_DELTABLAS}
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_blast}}, \code{\link{advanced_blast}}
#' @export
 

delta.blast <- function(query_file, 
                        subject_file,
                        cdd.path        = NULL,
                        seq_type        = "cds", 
                        format          = "fasta",
                        eval            = "1E-5",
                        max.target.seqs = 500,
                        path            = NULL, 
                        blast_params    = NULL,
                        comp_cores      = 1,
                        makedb_type     = "protein"){
        
        
        if(is.null(cdd.path))
                stop("Please specify the path to your DELTA-BLAST cdd database.")
        
       delta.output <- data.table::copy(
               
                        advanced_blast(
                            query_file      = query_file, 
                            subject_file    = subject_file, 
                            seq_type        = seq_type, 
                            format          = format,
                            path            = path, 
                            blast_params    = paste0("-evalue ",eval," -num_threads ",
                                                  comp_cores," -max_target_seqs ",
                                                  max.target.seqs," ",blast_params),
                            makedb_type     = makedb_type,
                            blast_algorithm = "deltablast",
                            db_path         = cdd.path )
                       
                       )
        
        return(delta.output)
        
}







