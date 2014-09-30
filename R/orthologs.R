#' @title Interface function to perform orthology inference for a set of organisms
#' @description This function takes nucleotide or protein sequences for a set of organisms 
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Further methods are: "BH" (BLAST best hit), "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), "OrthoMCL, "IP" (InParanoid).
#' @param path a character string specifying the path to the corresponding orthology inference tool.
#' For "BH" and "RBH": path to BLAST, "PO": path to ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL,
#' "IP": path to InParanoid.
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore computations.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @details This function takes sequence files of a query organism and a subject organism and performs orthology inference
#' using a defined orthology inference method to dectect orthologous genes.
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of orthologous genes
#' in the second column and the amino acid sequences in the third column.
#' @examples \dontrun{
#' 
#' ### Best Hit
#' 
#' # perform orthology inference using BLAST reciprocal best hit
#' # and fasta sequence files storing protein sequences
#' orthologs(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type = "protein", ortho_detection = "BH")
#'           
#' # multicore version          
#' orthologs(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type = "protein", ortho_detection = "BH", comp_cores = 2)          
#'           
#' # the same using CDS sequences
#' orthologs(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           seq_type = "cds", ortho_detection = "BH")
#' 
#' 
#' ### Reciprocal Best Hit
#' 
#' # perform orthology inference using BLAST reciprocal best hit
#' # and fasta sequence files storing protein sequences
#' orthologs(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type = "protein", ortho_detection = "RBH")
#'           
#' # multicore version          
#' orthologs(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type = "protein", ortho_detection = "RBH", comp_cores = 2)          
#'           
#' 
#' }
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}
#' @export
orthologs <- function(query_file,subject_file, seq_type = "protein",
                      format = "fasta",ortho_detection = "RBH", 
                      path = NULL, comp_cores = 1, quiet = FALSE){
        
        if(!is.element(ortho_detection, c("BH","RBH","PO","OrthoMCL","IP")))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        
        if(ortho_detection == "BH"){
                
                
                ortho_tbl <- data.table::copy(
                        blast_best(query_file = query_file, subject_file = subject_file, 
                                   path = path, comp_cores = comp_cores,
                                   seq_type = seq_type, format = format))
                
        }
        
        if(ortho_detection == "RBH"){
                
                
                ortho_tbl <- data.table::copy(
                        blast_rec(query_file = query_file, subject_file = subject_file, 
                                  path = path, comp_cores = comp_cores,
                                  seq_type = seq_type, format = format))
                
                
        }
        
        return(ortho_tbl)
        
}