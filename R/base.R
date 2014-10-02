#' @title Function to delete the internal folder hierarchy
#' @description This function deletes all internal folders that have been created
#' during pipeline processing. Internally this function uses \code{\link{unlink}}
#' to delete all folders created by the pipline.
#' @author Hajk-Georg Drost
#' @details Following folders are being removed: "_alignment", "_blast_db", and "_calculation".
#' @return Nothing is being returned.
#' @examples \dontrun{
#' 
#' # in case internal folders exist, they are being removed
#' clean_all_folders() 
#' 
#' }
#' @seealso \code{\link{divergence_stratigraphy}}
#' @export
clean_all_folders <- function(){
        
        if(file.exists("_alignment")){
                unlink("_alignment",recursive = TRUE, force = TRUE)
        }
        
        if(file.exists("_blast_db")){
                unlink("_blast_db",recursive = TRUE, force = TRUE)
        }
                
        
        if(file.exists("_calculation")){
                unlink("_calculation",recursive = TRUE, force = TRUE)
        }
                
}


set_path <- function(file, add.folder = NULL){
        
        f_sep <- .Platform$file.sep
        file_name <- unlist(strsplit(file,f_sep))
        file_name <- file_name[length(file_name)]
        curr_wd <- unlist(strsplit(getwd(),f_sep))
        wdir <- grepl(" ",curr_wd)
        curr_wd[wdir] <- stringr::str_replace(string = curr_wd[wdir],
                                              replacement = paste0("'",curr_wd[wdir],"'"), 
                                              pattern = curr_wd[wdir])
        
        curr_wd <- paste0(curr_wd,collapse = f_sep)
        
        if(!is.null(add.folder)){
                
                add.folder <- paste0(f_sep,add.folder)
                return(paste0(curr_wd,add.folder,file_name))
                
        } else {
                
                return(paste0(curr_wd,file_name))
        }
        
}


#' @title Function to filter the dNdS output of dNdS()
#' @description This function takes the output data.table returned by \code{\link{dNdS}} and
#' filters the output by the following criteria:
#' 
#' 1) all dN values having an NA value or omitted
#' 
#' 2) all dS values having an NA value or omitted
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
#' filter_dNdS( dNdS(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                   subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                   ortho_detection = "RBH", aa_aln_type = "multiple",
#'                   aa_aln_tool = "clustalw", aa_aln_path = "/path/to/clustalw/",
#'                   codon_aln_tool = "pal2nal", dnds_est.method = "YN", comp_cores = 1) , dnds.threshold = 2)
#' 
#' 
#' }
#' @import dplyr
#' @seealso \code{\link{divergence_stratum}}
#' @export
filter_dNdS <- function(dNdS_tbl,dnds.threshold = 2){
        
        return ( dNdS_tbl %>%
                        tbl_dt() %>%
                                dplyr::filter(!is.na(dN), !is.na(dS), dNdS <= dnds.threshold) )
        
}

