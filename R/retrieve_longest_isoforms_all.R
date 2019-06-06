#' @title Retrieve the longest isoforms of several proteome files stored in a folder
#' @description Based on a fasta file storing the peptide isoforms of gene loci and
#' an annotation file in \code{gtf} file format, this function extracts the longest
#' isoform per gene locus and stores the results in a new \code{fasta} file.
#' This procedure enables easier downstream analyses such as orthology inference etc
#' when dealing with proteome \code{fasta} files which usually include isoform peptides.
#' @param proteome_folder file path to proteome in \code{fasta} file format.
#' @param annotation_folder file path to the corresponding annotation file in \code{gtf} file format.
#' @param output_folder file path to new file storing only peptide sequences of the longest isoforms.
#' @param annotation_format format of \code{annotation_file}. Options are:
#' \itemize{
#' \item \code{annotation_file = "gff"} (default)
#' \item \code{annotation_file = "gtf"}
#' }
#' @author Hajk-Georg Drost
#' @examples \dontrun{
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
#' }
#' @export
retrieve_longest_isoforms_all <-
        function(proteome_folder,
                 annotation_folder,
                 output_folder,
                 annotation_format = "gff") {
                
                if (!file.exists(proteome_folder))
                        stop("Folder ", proteome_folder, " seems not to exist.", call. = FALSE)
                if (!file.exists(proteome_folder))
                        stop("Folder ", annotation_folder, " seems not to exist.", call. = FALSE)
                
                message("Starting retrieval of longest splice variants ...")
                message("Creating output folder '", output_folder,"' ...")
                dir.create(output_folder)
                proteome_list <- file.path(proteome_folder, list.files(proteome_folder))
                proteome_list <- proteome_list[!stringr::str_detect(proteome_list, "documentation")]
                annotation_list <- file.path(annotation_folder, list.files(annotation_folder))
                annotation_list <- annotation_list[!stringr::str_detect(annotation_list, "documentation")]
                
                for (i in seq_len(length(proteome_list))) {
                        org_name <- unlist(stringr::str_split(basename(proteome_list[i]), "[.]"))[1]
                        if (org_name != unlist(stringr::str_split(basename(annotation_list[i]), "[.]"))[1])
                                stop("The proteome file name ",basename(proteome_list[i])," and the annotation file name ", basename(annotation_list[1]), " seem not to match!",
                                     "Please make sure that organism file names *.faa and *.gtf/*.gff match so that the correct annotation file is matched with the correct proteome file.", call. = FALSE)
                        
                        message("Processing ", org_name, " ...")
                        retrieve_longest_isoforms(proteome_file = proteome_list[i], 
                                                  annotation_file = annotation_list[i], 
                                                  new_file = file.path(output_folder, paste0(org_name,".faa")))
                        
                }
                
                
                
                
                
        }