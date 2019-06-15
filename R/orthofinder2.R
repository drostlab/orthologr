#' @title Interface function to Orthofinder2
#' @description Run Orthofinder2 from R.
#' @param proteome_folder file path to a folder storing the proteome sequences of the species for which orthology inference shall be performed.
#' @param import_type type of \code{Orthofinder2} output that shall be imported after running \code{Orthofinder2}.
#' Options are:
#' \itemize{
#' \item \code{import_type = ""}
#' \item \code{import_type = ""}
#' \item \code{import_type = ""}
#' \item \code{import_type = ""}
#' }
#' @param comp_cores number of cores that shall be used for parallel processing. Default is \code{cores = 1}.
#' @param use_existing_output a logical value indicating whether or not an existing \code{Orthofinder2} 
#' output folder shall be used fo further import and processing. If \code{use_existing_output = TRUE}
#' is selected then please specify the file path to the \code{Orthofinder2} output folder in \code{proteome_folder}.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # specify species names
#' orgs <- c("Arabidopsis lyrata", 
#'           "Capsella rubella", "Solanum lycopersicum")
#' # download proteome files for all species          
#' biomartr::getProteomeSet(db = "refseq", organisms = orgs, path = "of_proteomes")
#' # download annotation files for all species          
#' biomartr::getGFFSet(db = "refseq", organisms = orgs, path = "of_gff")
#' # select longest splice variant per gene locus
#' retrieve_longest_isoforms_all(proteome_folder = "of_proteomes", 
#'                               annotation_folder = "of_gff",
#'                               annotation_format = "gff", 
#'                               output_folder = "of_proteomes_longest_sv")
#' # run orthofinder2 to infer ortho groups for the specified species
#' orthofinder2(proteome_folder = "of_proteomes_longest_sv", comp_cores = 4)
#' } 
#' @export
orthofinder2 <- function(proteome_folder, use_existing_output = FALSE, import_type = "", comp_cores = 1) {
        
        is_installed_orthofinder()
        
        ### add here a test that input sequences are amino acid sequences
        
        if (!file.exists(proteome_folder))
                stop("Please provide a valid path to the proteome folder. Your specified path '", proteome_folder, "' seems not to exist.", call. = FALSE)
        
        files <- list.files(proteome_folder)
        
        if (any(stringr::str_detect(files, "documentation"))) {
                message("Your species folder '", proteome_folder, " still contains the 'documentation' folder which will be moved from ", proteome_folder, " to ", getwd(), " to enable the Orthofinder2 search.")
                file.rename(file.path(proteome_folder, "documentation"), file.path(getwd(), "documentation"))
                unlink(file.path(proteome_folder, "documentation"), recursive = TRUE, force = TRUE)
        }
        
        if (!use_existing_output) {
                # determine the number of cores on a multicore machineYes
                cores <- parallel::detectCores()
                
                # in case one tries to use more cores than are available
                if (comp_cores > cores)
                        stop("You chose more cores than are available on your machine.", call. = FALSE)
                
                message("Running ", system("orthofinder", intern = TRUE)[1], " using ", comp_cores, " cores ...")
                message("Do your input fasta files store unique sequences per gene locus, e.g. the longest splice variant? If not, then have a look at ?orthologr::retrieve_longest_isoforms().")
                
                # output is stored in OrthoFinder/Results_DATE , where DATE is in format: Nov08 for 08 th November
                # or in this case OrthoFinder/Results_DATE_basename(proteome_folder)
                
                if (dirname(proteome_folder) == ".") {
                        system(paste0("orthofinder -f ", ws_wrap(file.path(getwd(), proteome_folder))," -t ", cores," -a ", cores," -S diamond -n ", basename(proteome_folder)))
                } else {
                        system(paste0("orthofinder -f ", ws_wrap(proteome_folder)," -t ", cores," -a ", cores," -S diamond -n ", basename(proteome_folder)))
                }
                
                message("Orthofinder2 finished successfully and stored all results in ", ifelse(dirname(proteome_folder) == ".", ws_wrap(file.path(getwd(), proteome_folder)), ws_wrap(proteome_folder)))
                
                if (import_type == "") {
                        
                }
                
                
        } else {
                
                
        }
        
        
        
}
