#' @title Filter dNdS values
#' @description This function takes the output data.table returned by \code{\link{dNdS}} and
#' filters the output by the following criteria:
#' 
#' 1) all dN values having an NA value are omitted
#' 
#' 2) all dS values having an NA value are omitted
#' 
#' 3) all dNdS values >= the specified \code{dnds.threshold} are omitted
#' 
#' @param dNdS_tbl a \code{data.table} returned by \code{\link{dNdS}}.
#' @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @details
#' 
#' The dNdS ratio quantifies the selection pressure acting on a given protein sequence.
#' 
#' It is proposed that:
#' 
#'  1) dNdS values < 1 reflect stabilizing selection of a protein
#'  
#'  2) dNdS value = 1 reflect neutral selection of a protein
#'  
#'  3) dNdS values > 1 reflect variational selection of a protein
#' 
#' Now an assumption must be generated to allow for dNdS filtering.
#' Given a \code{dnds.threshold} all value above this threshold are omitted from the dataset.
#' 
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' 
#' filter_dNdS( dNdS( query_file = 
#'                    system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                    subject_file = 
#'                    system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                    ortho_detection = "RBH", 
#'                    aa_aln_type = "pairwise",
#'                    aa_aln_tool = "NW",
#'                    codon_aln_tool = "pal2nal", 
#'                    dnds_est.method = "Comeron", 
#'                    comp_cores = 1), 
#'              dnds.threshold = 2)
#' 
#' # a small example using clustalw
#' filter_dNdS( dNdS( query_file      = 
#'                    system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                    subject_file    = 
#'                    system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                    ortho_detection = "RBH", 
#'                    aa_aln_type     = "pairwise",
#'                    aa_aln_tool     = "NW",
#'                    codon_aln_tool  = "pal2nal", 
#'                    dnds_est.method = "Comeron", 
#'                    comp_cores      = 1) , 
#'                dnds.threshold  = 2)
#'                   
#' }
#' @seealso \code{\link{divergence_stratigraphy}}
#' @export

filter_dNdS <- function(dNdS_tbl,dnds.threshold = 2){
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        dN <- dS <- NULL
        
        return( dplyr::filter(dNdS_tbl,!is.na(dN), !is.na(dS), dNdS <= dnds.threshold) )
        
}