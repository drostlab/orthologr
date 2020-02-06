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
#' \item \code{annotation_file = "gff"} (default)
#' \item \code{annotation_file = "gtf"} 
#' }
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # retrieve example data from ENSEMBLGENOMES
#' proteome <- biomartr::getProteome(db = "refseq", organism = "Arabidopsis thaliana")
#' annotation <- biomartr::getGFF(db = "refseq", organism = "Arabidopsis thaliana")
#' # retrieve longest isoforms and store in new file
#' retrieve_longest_isoforms(proteome_file = proteome, 
#'                           annotation_file = annotation, 
#'                           new_file = "Athaliana_pep_longest.fa")
#' # import new file into R session                          
#' Athaliana_pep_longest <- Biostrings::readAAStringSet("Athaliana_pep_longest.fa")
#' }
#' @export
retrieve_longest_isoforms <- function(proteome_file, annotation_file, new_file, annotation_format = "gff") {
        
        if (!file.exists(proteome_file))
                stop("File ", proteome_file, " seems not to exist.", call. = FALSE)
        
        if (!file.exists(annotation_file))
                stop("File ", annotation_file, " seems not to exist.", call. = FALSE)
        
        if (!is.element(annotation_format, c("gtf", "gff")))
                stop("Please specify a valid annotation format. Options are: annotation_file = 'gtf' or annotation_file = 'gff'.", call. = FALSE)
        
        '.' <- NULL
        
        message("Extracting longest isoform from ", basename(proteome_file), " ...")
        message("Importing proteome file ", basename(proteome_file), " ...")
        pep_file <- Biostrings::readAAStringSet(proteome_file)
        new_names <- unlist(lapply(pep_file@ranges@NAMES, function(x) {
                unlist(stringr::str_split(x," "))[1]
        }))
        names(pep_file) <- new_names
        
        # check if annotation file is corrupt and remove outliers
        annotation_file <- check_annotation(annotation_file, remove_annotation_outliers = TRUE)
        
        if (annotation_format == "gtf") {
                
                pep_file_tibble <- tibble::tibble(transcript_id = new_names, width = pep_file@ranges@width)
                #print(pep_file_tibble)
                
                message("Importing gtf file ", basename(annotation_file), " ...")
                gtf_import <- tibble::as_tibble(rtracklayer::import(annotation_file))
                message("Filter for gene_biotype == 'protein_coding' AND type == 'transcript' ...")
                gene_biotype <- gene_id <- transcript_id <- type <- width <- NULL
                gtf_import <- dplyr::filter(gtf_import, gene_biotype == "protein_coding", type == "transcript")
                message("After filtering, the gtf file includes ", nrow(gtf_import), " rows.")
                gtf_import <- dplyr::select(gtf_import, gene_id, transcript_id)
                #print(gtf_import)
                # join transcript_ids from annotation file with peptide sequence file
                message("Join transcript_ids from FASTA header with transcript_ids from GTF annotation file ...")
                pep_file_tibble_joined <- dplyr::inner_join(pep_file_tibble, gtf_import, by = "transcript_id")
                message("The joined file contains ", nrow(pep_file_tibble_joined), " rows.")
                message("Removing rows with NA's in columns 'transcript_id', 'width', and 'gene_id'")
                pep_file_tibble_joined <- dplyr::filter(pep_file_tibble_joined, !is.na(transcript_id), !is.na(width), !is.na(gene_id))
                message("After NA removal the file contains ", nrow(pep_file_tibble_joined), " rows.")
                
                #print(pep_file_tibble_joined)
                extract_longest_transcript <- function(x) {
                        pos <- which.max(x$width)[1]
                        return(dplyr::slice(x, pos))
                }
                message("Retrieving longest isoforms ...")
                res <- dplyr::do(dplyr::group_by(pep_file_tibble_joined, gene_id), extract_longest_transcript(.))
                matched_ids <- stats::na.omit(match(res$transcript_id, pep_file@ranges@NAMES))
                
                if (length(matched_ids) == 0)
                        stop("Gene_ID or Peptide_IDs between the gtf file and the headers of the fasta file did not match!",
                             " Please make sure that the 'transcript_id' column in the gtf file matches the IDs given in the headers of the fasta file.", call. = FALSE)
                if (length(matched_ids) != length(unique(matched_ids)))
                        stop("There should only be unique gene_ids in the file, but non-unique IDs were found! Please check what could have gone wrong. with the input files.", call. = FALSE)
                
                message("\n")
                message("Writing ", length(matched_ids), " unique longest peptide sequences (from initially ", length(pep_file), " isoforms) to new fasta file ", new_file, " ...")
                Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
                message("Retrieval finished successfully.")
        }
        
        if (annotation_format == "gff") {
                
                pep_file_tibble <- tibble::tibble(protein_id = new_names, width = pep_file@ranges@width)
                #print(pep_file_tibble)
                message("Importing gff file ", basename(annotation_file), " ...")
                gff_import <- tibble::as_tibble(rtracklayer::import(annotation_file))
                gene_biotype <- gene <- transcript_id <- type <- width <- NULL
                gff_import <- dplyr::filter(gff_import, type == "CDS")
                protein_id <- NULL
                gff_import <- dplyr::select(gff_import, gene, protein_id)
                #print(gff_import)
                # join transcript_ids from annotation file with peptide sequence file
                pep_file_tibble_joined <- dplyr::inner_join(pep_file_tibble, gff_import, by = "protein_id")   
                pep_file_tibble_joined <- dplyr::filter(pep_file_tibble_joined, !is.na(protein_id), !is.na(width))
                #print(pep_file_tibble_joined)
                extract_longest_transcript <- function(x) {
                        pos <- which.max(x$width)[1]
                        return(dplyr::slice(x, pos))
                }
                message("Retrieving longest isoforms ...")
                res <- dplyr::do(dplyr::group_by(pep_file_tibble_joined, gene), extract_longest_transcript(.))
                matched_ids <- stats::na.omit(match(res$protein_id, pep_file@ranges@NAMES))
                
                if (length(matched_ids) == 0)
                        stop("Gene_ID or Peptide_IDs between the gtf file and the headers of the fasta file did not match!",
                             " Please make sure that the 'protein_id' column in the gtf file matches the IDs given in the headers of the fasta file.", call. = FALSE)
                
                if (length(matched_ids) != length(unique(matched_ids)))
                        stop("There should only be unique gene_ids in the file, but non-unique IDs were found! Please check what could have gone wrong. with the input files.", call. = FALSE)
                
                message("\n")
                message("Writing ", length(matched_ids), " unique longest peptide sequences (from initially ", length(pep_file), " isoforms) to new fasta file ", new_file, " ...")
                Biostrings::writeXStringSet(pep_file[matched_ids], new_file)
                message("Retrieval finished successfully.")
        }
        
        
}
