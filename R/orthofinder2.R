#' @title Interface function to Orthofinder2
#' @description Run Orthofinder2 from R.
#' @param proteome_folder file path to a folder storing the proteome sequences of the species for which orthology inference shall be performed.
#' @param comp_cores number of cores that shall be used for parallel processing. Default is \code{cores = 1}.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # specify species names
#' orgs <- c("Arabidopsis thaliana", "Arabidopsis lyrata", 
#'           "Capsella rubella", "Solanum lycopersicum")
#' # download proteome files for all species          
#' biomartr::getProteomeSet(db = "refseq", organisms = orgs, path = "of_proteomes")
#' # download annotation files for all species          
#' biomartr::getGFFSet(db = "refseq", organisms = orgs, path = "of_gff")
#' # select longest splice variant per gene locus
#' dir.create("of_proteomes_longest_sv")
#' retrieve_longest_isoforms(proteome_file = "of_proteomes/Alyrata.faa", 
#'                           annotation_file = "of_gff/Alyrata.gff",
#'                           annotation_format = "gff", 
#'                           new_file = "of_proteomes_longest_sv/Alyrata.faa")
#' # run orthofinder2 to infer ortho groups for the specified species
#' orthofinder2(proteome_folder = "of_proteomes", comp_cores = 4)
#' } 
#' @export
orthofinder2 <- function(proteome_folder, comp_cores = 1) {
        
        is_installed_orthofinder()
        
        ### add here a test that input sequences are amino acid sequences
        
        if (!file.exists(proteome_folder))
                stop("Please provide a valid path to the proteome folder. Your specified path '", proteome_folder, "' seems not to exist.", call. = FALSE)
        
        files <- list.files(proteome_folder)
        
        if (stringr::str_detect(files, "documentation")) {
                message("Your species folder '", proteome_folder, " still contains the 'documentation' folder which will be moved from ", proteome_folder, " to ", getwd(), " to enable the Orthofinder2 search.")
                file.rename(file.path(proteome_folder, "documentation"), file.path(getwd(), "documentation"))
                unlink(file.path(proteome_folder, "documentation"), recursive = TRUE, force = TRUE)
        }
        
        # determine the number of cores on a multicore machineYes
        cores <- parallel::detectCores()
        
        # in case one tries to use more cores than are available
        if (comp_cores > cores)
                stop("You chose more cores than are available on your machine.", call. = FALSE)
        
        message("Running ", system("orthofinder", intern = TRUE)[1], " using ", comp_cores, " cores ...")
        message("Do your input fasta files store unique sequences per gene locus, e.g. the longest splice variant? If not, then have a look at ?orthologr::retrieve_longest_isoforms().")
        
        # output is stored in OrthoFinder/Results_DATE , where DATE is in format: Nov08 for 08 th November
        # or in this case OrthoFinder/Results_DATE_basename(proteome_folder)
        system(paste0("orthofinder -f ", proteome_folder," -t ", cores," -a ", cores," -S diamond -n ", basename(proteome_folder)))
}