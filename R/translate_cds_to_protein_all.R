#' @title Translate coding sequences into amino acid sequences for multiple files
#' @description A helper function that takes \code{fasta} files storing coding sequences
#' as input and translates these coding sequences into amino acid sequences
#' storing them as \code{fasta} output files.
#' @param input_folder file path to folder storing the coding sequences \code{fasta} files.
#' @param output_folder name or file path to a folder that that shall be generated to store the output \code{fasta} files.
#' @param delete_corrupt_cds delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @author Hajk-Georg Drost
#' @examples dontrun{
#' # install.packages("biomartr")
#' # download coding sequences of Arabidopsis thaliana, Arabidopsis lyrata, and Capsella rubella
#' org_list <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella")
#' biomartr::getCDSSet(db = "refseq",
#'              organism = org_list,
#'              gunzip = TRUE,
#'              path = "cds_examples")
#' # translate coding sequences into amino acid sequences
#' translate_cds_to_protein_all(input_folder = "cds_examples", 
#'                              output_folder = "translated_seqs",
#'                              delete_corrupt_cds = FALSE)
#' }
#' @export
translate_cds_to_protein_all <-
        function(input_folder,
                 output_folder,
                 delete_corrupt_cds = FALSE) {
                
                if (!file.exists(input_folder))
                        stop(
                                "The folder '",
                                input_folder,
                                "' seems not to exist. Please check your file path.",
                                call. = FALSE
                        )
                
                if (!file.exists(output_folder)) {
                        message("Creating new output folder '", output_folder, "'.")
                        dir.create(output_folder)
                }
                
                file_list <- list.files(input_folder)
                
                found_documentation <- stringr::str_detect(file_list, "documentation")
                if (length(found_documentation) > 0)
                        file_list <- file_list[!found_documentation]
                
                if (length(file_list) == 0)
                        stop("No files were found in your specified input_folder.", call. = FALSE)
                
                for (i in seq_len(length(file_list))) {
                        message("Processing file ", file_list[i], " ...")
                        translate_cds_to_protein(
                                input_file = file.path(input_folder, file_list[i]),
                                output_file = file.path(output_folder, file_list[i]),
                                delete_corrupt_cds = delete_corrupt_cds
                        )
                        message("\n")
                }
                        
                message("All input files were successfully translated and stored at ", output_folder,".")
        }