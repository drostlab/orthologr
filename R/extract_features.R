#' @title Helper function to extract gene loci and splice variant IDs from GFF files
#' @description Extract gene loci and splice variant IDs from GFF files.
#' @param x a data.frame converted from a \code{gff} or \code{gtf} file.
#' @param format either \code{format = "gtf"} or \code{format = "gff"}.
#' @export
extract_features <- function(x, format) {
        
        if (!is.element(format, c("gtf", "gff")))
                stop("Please choose a format that is supported by this function: format = 'gtf' or format = 'gff'.", call. = FALSE)
        
        if (format == "gtf") {
                if (!all(c("gene", "CDS") %in% names(x)))
                        stop("Please make sure that the gene and CDS biotypes in your gtf file are labeled as 'gene' and 'CDS'.", call. = FALSE)
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "CDS")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$gene_id, length(unique(get_mRNA_id$transcript_id))),
                                              query_id = unique(get_mRNA_id$transcript_id))
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        if (format == "gff") {
                if (!all(c("gene", "mRNA") %in% names(x)))
                        stop("Please make sure that the gene and mRNA biotypes in your gtf file are labeled as 'gene' and 'mRNA'.", call. = FALSE)
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "mRNA")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$ID, length(unique(get_mRNA_id$ID))),
                                              query_id = unique(get_mRNA_id$ID))
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        gene_locus_id <- query_id <- NULL
        res <- dplyr::filter(res, !is.na(gene_locus_id), !is.na(query_id))
        return(res)
}