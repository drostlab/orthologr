#' @title Helper fucntion to extract gene locu and splice variant IDs from GFF files
extract_features <- function(x, format) {
        
        if (!is.element(format, c("gtf", "gff")))
                stop("Please choose a format that is supported by this function: format = 'gtf' or format = 'gff'.", call. = FALSE)
        
        if (format == "gtf") {
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "CDS")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$gene_id, length(get_mRNA_id$transcript_id)),
                                              query_id = get_mRNA_id$transcript_id)
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        if (format == "gff") {
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "mRNA")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$ID, length(get_mRNA_id$ID)),
                                              query_id = get_mRNA_id$ID)
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        
        if (format == "gff") {
                get_gene_id <- dplyr::filter(x, type == "gene")
                get_mRNA_id <- dplyr::filter(x, type == "CDS")
                
                if (nrow(get_gene_id) > 0 & nrow(get_mRNA_id) > 0) {
                        
                        res <- tibble::tibble(gene_locus_id = rep(get_gene_id$ID, length(get_mRNA_id$ID)),
                                              query_id = get_mRNA_id$ID)
                        
                } else {
                        res <- tibble::tibble(gene_locus_id = NA, query_id = NA)
                }
        }
        return(res)
}