#' @title Function to perform divergence stratigraphy
#' @description This function takes a query organism and performs
#' divergence stratigraphy (Quint et al.,2012 ; Drost et al. 2014) against a
#' closely related subject organism.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Available methods are: "BH" (BLAST best hit), "RBH" (BLAST reciprocal best hit).
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param mafft_path a character string specifying the path to the multiple alignment program MAFFT (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform
#'  parallel computations on a multicore machine.
#'  @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out to the console.
#' @param clean_folders a logical value specifying whether the internal folder structure shall be deleted (cleaned) after
#'  processing this function. Default is \code{clean_folders} = \code{FALSE}.
#' @details Introduced by Quint et al.,2012 and extended in Drost et al. 2014, divergence stratigraphy
#'  is the process of quantifying the selection pressure (in terms of amino acid sequence divergence) acting on
#'  orthologous genes between closely related species. The resulting sequence divergence map (short divergence map),
#'  stores the divergence stratum in the first column and the query_id of inferred orthologous genes in the second column.
#'  
#'  Following steps are performed to obtain a standard divergence map based on divergence_stratigraphy:
#'  
#'  1) Orthology Inference using BLAST best hit ("BH") or BLAST reciprocal best hit ("RBH")
#'  
#'  2) Pairwise amino acid alignments of orthologous genes using the MAFFT program
#'  
#'  3) Codon alignments of orthologous genes using PAL2NAL
#'  
#'  4) dNdS estimation using Comeron's method (1995)
#'  
#'  5) Assigning dNdS values to divergence strata (deciles)
#'  
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @references
#'  
#'  Quint M et al. (2012). "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#'  
#'  Drost HG et al. (2014). "Active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis".
#'  
#' @examples \dontrun{
#'  
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                          eval = "1E-5", ortho_detection = "BH",mafft_path = "path/to/mafft",
#'                          comp_cores = 1, quiet = TRUE, clean_folders = TRUE)
#'  
#'  
#'  }
#'  
#' @return A data.table storing the divergence map of the query organism.
#' @seealso \code{\link{dNdS}}, \code{\link{substitutionrate}}, \code{\link{multi_aln}},
#'   \code{\link{codon_aln}}, \code{\link{blast_best}}, \code{\link{blast_rec}}
#' @export
divergence_stratigraphy <- function(query_file, subject_file, eval = "1E-5",
                                    ortho_detection = "BH", blast_path = NULL, 
                                    mafft_path = NULL, comp_cores = 1,dnds.threshold = 2,
                                    quiet = FALSE, clean_folders = FALSE){
        
        if(!is.ortho_detection_method(ortho_detection))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        dNdS_tbl <- filter_dNdS( dNdS(query_file = query_file,
                                      subject_file = subject_file,
                                      ortho_detection = ortho_detection,
                                      aa_aln_type = "multiple", aa_aln_tool = "mafft",
                                      aa_aln_path = mafft_path, codon_aln_tool = "pal2nal", 
                                      dnds_est.method = "Comeron", comp_cores = comp_cores, 
                                      quiet = quiet) , dnds.threshold = dnds.threshold)
        
        # divergence map: standard = col1: divergence stratum, col2: query_id
        dm_tbl <- DivergenceMap( dNdS_tbl[ ,list(dNdS,query_id)] )
        
        
        if(clean_folders)
                clean_all_folders()
        
        return ( dm_tbl )
        
}



#' @title Function to sort dNdS values into divergence strata.
#' @description This function takes a data.table returned by dNdS
#' and sorts the corresponding dNdS value into divergence strata (deciles).
#' @param dNdS_tbl a data.table object returned by \code{\link{dNdS}}, but only selecting dNdS_obj[ , list(dNdS,query_id)].
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @examples \dontrun{
#' 
#' # get a divergence map of example sequences
#' dNdS_tbl <-  dNdS( query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                    subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                    ortho_detection = "RBH", aa_aln_type = "multiple",
#'                    aa_aln_tool = "clustalw", codon_aln_tool = "pal2nal", 
#'                    dnds_est.method = "YN", comp_cores = 1 )
#'                  
#' divMap <- DivergenceMap( dNdS_tbl[ , list(dNdS,query_id)] )                       
#' 
#'               
#' 
#' }
#' @seealso \code{\link{divergence_stratigraphy}}
#' @return a data.table storing a standard divergence map.
#' @export
DivergenceMap <- function(dNdS_tbl){
        
        DecileValues <- stats::quantile(dNdS_tbl[ , dNdS],probs = seq(0, 1, 0.1))
        
        for(i in length(DecileValues):2){
                
                AllGenesOfDecile_i <- na.omit(which((dNdS_tbl[ , dNdS] < DecileValues[i]) & (dNdS_tbl[ , dNdS] >= DecileValues[i-1])))
                dNdS_tbl[AllGenesOfDecile_i, dNdS:=(i-1)] 
                
        }
        
        ## assigning all KaKs values to Decile-Class : 10 which have the exact Kaks-value
        ## as the 100% quantile, because in the loop we tested for < X% leaving out
        ## the exact 100% quantile
        dNdS_tbl[which(dNdS_tbl[ , dNdS] == DecileValues[length(DecileValues)]) , 1] <- 10
        
        return(dNdS_tbl)
        
}





