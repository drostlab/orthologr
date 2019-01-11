#' @title Retrieve the longest isoforms from a proteome file and save results as fasta file
#' @description Based on a fasta file storing the peptide isoforms of gene loci and
#' an annotation file in \code{gtf} file format, this function extracts the longest
#' isoform per gene locus and stores the results in a new \code{fasta} file.
#' This procedure enables easier downstream analyses such as orthology inference etc
#' when dealing with proteome \code{fasta} files which usually include isoform peptides.
#' @param proteome_file file path to proteome in \code{fasta} file format.
#' @param annotation_file file path to the corresponding annotation file in \code{gtf} file format.
#' @param new_file file path to new file storing only peptide sequences of the longest isoforms.
#' @param annotation_format format of \code{annotation_file}. Options are:
#' \itemize{
#' \item \code{annotation_file = "gtf"} (default)
#' \item \code{annotation_file = "gff"}
#' }
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # retrieve example data from ENSEMBLGENOMES
#' proteome <- biomartr::getProteome(db = "ensemblgenomes", organism = "Arabidopsis thaliana")
#' biomartr::getGTF(db = "ensemblgenomes", organism = "Arabidopsis thaliana")
#' R.utils::gunzip("ensembl/annotation/Arabidopsis_thaliana.TAIR10.42_ensemblgenomes.gtf.gz")
#' annotation <- "ensembl/annotation/Arabidopsis_thaliana.TAIR10.42_ensemblgenomes.gtf"
#' # retrieve longest isoforms and store in new file
#' retrieve_longest_isoforms(proteome_file = proteome, 
#'                           annotation_file = annotation, 
#'                           new_file = "Athaliana_pep_longest.fa")
#' # import new file into R session                          
#' Athaliana_pep_longest <- Biostrings::readAAStringSet("Athaliana_pep_longest.fa")
#' }
#' @export
retrieve_longest_isoforms <- function(proteome_file, annotation_file, new_file, annotation_format = "gtf") {
        
        if (!fs::file_exists(proteome_file))
                stop("File ", proteome_file, " seems to not exist.", call. = FALSE)
        if (!fs::file_exists(annotation_file))
                stop("File ", annotation_file, " seems to not exist.", call. = FALSE)
        if (!is.element(annotation_file, c("gtf", "gff")))
                stop("Please specify a valid annotation_file. Options are: annotation_file = 'gtf' or annotation_file = 'gff'.", call. = FALSE)
        
        message("Extracting longest isoform from ", basename(proteome_file), " ...")
        message("Importing proteome file ", basename(proteome_file), " ...")
        pep_file <- Biostrings::readAAStringSet(proteome_file)
        new_names <- unlist(lapply(pep_file@ranges@NAMES, function(x) {
                unlist(stringr::str_split(x," "))[1]
        }))
        names(pep_file) <- new_names
        pep_file_tibble <- tibble::tibble(transcript_id = new_names, width = pep_file@ranges@width)
        #print(pep_file_tibble)
        
        if (annotation_file == "gtf") {
                message("Importing gtf file ", basename(annotation_file), " ...")
                gtf_import <- tibble::as_tibble(rtracklayer::import(annotation_file))
                gene_biotype <- gene_id <- transcript_id <- type <- width <- NULL
                gtf_import <- dplyr::filter(gtf_import, gene_biotype == "protein_coding", type == "transcript")
                gtf_import <- dplyr::select(gtf_import, gene_id, transcript_id)
                #print(gtf_import)
                # join transcript_ids from annotation file with peptide sequence file
                pep_file_tibble_joined <- dplyr::inner_join(pep_file_tibble, gtf_import, by = "transcript_id")   
                pep_file_tibble_joined <- dplyr::filter(pep_file_tibble_joined, !is.na(transcript_id), !is.na(width), !is.na(gene_id))
                #print(pep_file_tibble_joined)
                extract_longest_transcript <- function(x) {
                        pos <- which.max(x$width)[1]
                        return(dplyr::slice(x, pos))
                }
                message("Retrieving longest isoforms ...")
                res <- dplyr::do(dplyr::group_by(pep_file_tibble_joined, gene_id), extract_longest_transcript(.))
                matched_ids <- na.omit(match(res$transcript_id, pep_file@ranges@NAMES))
                
                if (length(matched_ids) == 0)
                        stop("Gene_ID or Peptide_IDs between the gtf file and the headers of the fasta file did not match!",
                             " Please make sure that the 'transcript_id' column in the gtf file matches the IDs given in the headers of the fasta file.", call. = FALSE)
                
                message("\n")
                message("Writing ", length(matched_ids), " longest peptide sequences (from initially ", length(pep_file), " isoforms) to new fasta file ", new_file, " ...")
                Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
                message("Retrieval finished successfully.")
        }
        
}
