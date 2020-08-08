#' @title Orthology Inference of lncRNAs
#' @description This function takes nucleotide or protein sequences for a set of organisms 
#' and performs orthology inference to detect orthologous genes within the given organisms
#' based on selected orthology inference programs.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the paths to the sequence files of interest (subject organisms).
#' @param task nucleotide search task option. Options are:
#' \itemize{
#' \item \code{task = "blastn"} : Standard nucleotide-nucleotide comparisons (default) - Traditional BLASTN requiring an exact match of 11.
#' \item \code{task = "blastn-short"} : Optimized nucleotide-nucleotide comparisons for query sequences shorter than 50 nucleotides.
#' \item \code{task = "dc-megablast"} : Discontiguous megablast used to find somewhat distant sequences.
#' \item \code{task = "megablast"} : Traditional megablast used to find very similar (e.g., intraspecies or closely related species) sequences.
#' \item \code{task = "rmblastn"}
#' }
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Options are:
#' \itemize{
#' \item \code{ortho_detection = "RBH"} (BLAST reciprocal best hit) (Default)
#' \item \code{ortho_detection = "BH"} (BLAST best hit)
#' }
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep. Please be aware that max.target.seqs selects best hits based on the database entry and not by the best e-value. 
#' See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param output.path path to which output shall be stored.
#' @param path a character string specifying the path to the corresponding orthology inference tool.
#' For "BH" and "RBH": path to BLAST, "PO": path to ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL.
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore computations.
#' @details 
#' This function takes sequence files of a query organism and a subject organism and performs orthology inference
#' using a defined orthology inference method to dectect orthologous genes.
#' 
#' The following interfaces are implemented in the \code{orthologs} function:
#'  
#' BLAST based methods:
#' 
#' \itemize{
#'   \item BLAST best hit (BH) 
#'   \item BLAST reciprocal best hit (RBH)
#' }
#' 
#' 
#' @author Hajk-Georg Drost
#' @return A data.table storing the query_ids of orthologous genes in the first column, the subject_ids of orthologous genes
#' in the second column and the amino acid sequences in the third column.
#' @references
#' 
#' BLAST: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
#' 
#' ProteinOrtho: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
#'
#' @examples \dontrun{
#' 
#' ### BLAST Reciprocal Best Hit
#' # perform orthology inference using BLAST reciprocal best hit
#' # and fasta sequence files storing protein sequences
#' orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           ortho_detection = "RBH")
#'           
#'           
#' ### BLAST Best Hit
#' # perform orthology inference using BLAST best hit
#' # and fasta sequence files storing protein sequences
#' orthologs.lnc(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           ortho_detection = "BH")
#' 
#'
#' # multicore version          
#' orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'           ortho_detection = "RBH", 
#'           comp_cores      = 2)          
#'           
#'           
#'           
#' }
#' @seealso \code{\link{blast_rec}}, \code{\link{dNdS}}
#' @export

orthologs_lnc <- function(query_file,
                      subject_file,
                      task            = "blastn",
                      eval            = "1E-5", 
                      ortho_detection = "RBH",
                      max.target.seqs = 10000,
                      output.path = getwd(),
                      comp_cores      = 1,
                      path = NULL) {
        
        message("Starting orthology inference of lncRNAs between query species: ", basename(query_file), " and subject species: ", basename(subject_file))
        
        if (!is.element(ortho_detection,
                       c("BH", "RBH")))
                stop("Please choose a orthology detection method that is supported by this function.", call. = FALSE)
        
        if (ortho_detection == "BH") {
                
                if (length(subject_file) > 1)
                        stop("The BLAST best hit method is only defined for pairwise comparisons.", call. = FALSE)
                
                ortho_tbl <- metablastr::blast_best_hit(
                        query      = query_file,
                        subject =  subject_file,
                        search_type = "nucleotide_to_nucleotide",
                        cores = comp_cores,
                        task = task,
                        evalue = eval, 
                        output.path = tempdir(),
                        max.target.seqs = max.target.seqs,
                        blast.path = path
                )
                
        }

        if (ortho_detection == "RBH") {
                
                if (length(subject_file) > 1)
                        stop("The BLAST best reciprocal hit method is only defined for pairwise comparisons.", call. = FALSE)
                
                ortho_tbl <- metablastr::blast_best_reciprocal_hit(
                        query      = query_file,
                        subject =  subject_file,
                        search_type = "nucleotide_to_nucleotide",
                        cores = comp_cores,
                        task = task,
                        evalue = eval, 
                        output.path = tempdir(),
                        max.target.seqs = max.target.seqs,
                        blast.path = path
                )

        }
        
        message("Orthology inference finished successsfully. BLAST output tables were stored in at tempdir() -> ", tempdir(), " .")
        return(ortho_tbl)
}