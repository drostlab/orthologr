#' @title Interface function to Orthofinder2
#' @description Run Orthofinder2 from R.
#' @param proteome_folder file path to a folder storing the proteome sequences of the species for which orthology inference shall be performed.
#' @param comp_cores number of cores that shall be used for parallel processing. Default is \code{cores = 1}.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # download proteome files
#' biomartr::getProteome(db = "refseq", organism = "Arabidopsis thaliana", path = "of_proteomes")
#' biomartr::getProteome(db = "refseq", organism = "Arabidopsis lyrata", path = "of_proteomes")
#' biomartr::getProteome(db = "refseq", organism = "Capsella rubella", path = "of_proteomes")
#' biomartr::getProteome(db = "refseq", organism = "Solanum lycopersicum", path = "of_proteomes")
#' 
#' } 
#' @export
orthofinder2 <- function(proteome_folder, comp_cores = 1) {
        
        is_installed_orthofinder()
        
        if (!file.exists(proteome_folder))
                stop("Please provide a valid path to the proteome folder. Your specified path '", proteome_folder, "' seems not to exist.", call. = FALSE)
        
        # determine the number of cores on a multicore machine
        cores <- parallel::detectCores()
        
        # in case one tries to use more cores than are available
        if (comp_cores > cores)
                stop("You chose more cores than are available on your machine.")
        
        message("Running ", system("orthofinder", intern = TRUE)[1])
        
        # output is stored in OrthoFinder/Results_DATE , where DATE is in format: Nov08 for 08 th November
        # or in this case OrthoFinder/Results_DATE_basename(proteome_folder)
        system(paste0("orthofinder -f ", proteome_folder," -t ", cores," -a ", cores," -S diamond -n ", basename(proteome_folder)))
}