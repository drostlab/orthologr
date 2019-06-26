#' @title Translate coding sequences into amino acid sequences
#' @description A helper function that takes a \code{fasta} file storing coding sequences
#' as input and translates these coding sequences into amino acid sequences
#' storing them as \code{fasta} output file.
#' @param input_file file path to \code{fasta} file storing coding sequences (DNA).
#' @param output_file name or file path in which translated amino acid sequences (AA)
#' shall be stored as \code{fasta} file.
#' @param delete_corrupt_cds delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @author Hajk-Georg Drost
#' @examples dontrun{
#' # install.packages("biomartr")
#' # download coding sequences of Arabidopsis thaliana
#' Ath_file_path <- biomartr::getCDS(db = "refseq",
#'                      organism = "Arabidopsis thaliana",
#'                      gunzip = TRUE)
#' # translate coding sequences into amino acid sequences
#' translate_cds_to_protein(Ath_file_path, "Ath_aa_seqs.fasta")
#' }
#' @export
translate_cds_to_protein <-
        function(input_file,
                 output_file,
                 delete_corrupt_cds = FALSE) {
                if (!file.exists(input_file))
                        stop(
                                "The file '",
                                input_file,
                                "' seems not to exist. Please check your file path.",
                                call. = FALSE
                        )
                message("Starting translation of file ",
                        basename(input_file),
                        " ...")
                cds_file <-
                        read.cds(
                                file = input_file,
                                format = "fasta",
                                delete_corrupt_cds = delete_corrupt_cds
                        )
                cds_seqs <- as.data.frame(cds_file)
                seq_vector <- cds_seqs$seqs
                names(seq_vector) <- cds_seqs$geneids
                seqs_biostrings <- Biostrings::DNAStringSet(seq_vector, use.names = TRUE)   
                protein_file <-
                        Biostrings::translate(seqs_biostrings, if.fuzzy.codon = "solve")
                Biostrings::writeXStringSet(protein_file, output_file)
                message("Storing translated file at ", output_file, ".")
        }