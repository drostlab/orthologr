#' @title Select orthologs based on either gene locus or splice variant
#' @description This function ...
#' @param dnds_tbl
#' @param annotation_file_query file path to query annotation file in either \code{gtf} or \code{gff} format.
#' @param annotation_file_subject file path to query annotation file in either \code{gtf} or \code{gff} format.
#' @param collapse_by locus type by which orthologs should be determined. Options are:
#' \itemize{
#' \item \code{type = "gene_locus"}: the splice variant having the lowest e-value of the BLAST search is chosen to represent the gene locus for the orthology relationship.
#' \item \code{type = "splice_variant"}: all splice variants of homologous gene loci will be returned.
#' }
#' @param format a vector of length 2 storing the annotation file formats of the query annotation file and subject annotation file: either \code{gtf} or \code{gff} format. E.g. \code{format = c("gtf","gtf")}.
#' @author Hajk-Georg Drost
#' @export
select_orthologs <-
        function(dnds_tbl,
                 annotation_file_query,
                 annotation_file_subject,
                 collapse_by = "gene_locus",
                 format = c("gtf","gtf")) {
                
                if (!is.element(collapse_by, c("gene_locus", "splice_variant")))
                        stop(
                                "Please specify a collapse method provided by this function: collapse_by = 'gene_locus' or collapse_by = 'splice_variant'.",
                                call. = FALSE
                        )
                
                if (!all(is.element(format, c("gtf", "gff"))))
                        stop("Please provide a vector of length 2 that contains available formats: either 'gtf' or 'gff'.", call. = FALSE)
                
                if (length(format) != 2)
                        stop("Please provide a vector of length 2 to specify the query and subject annotation file formats. E.g. format = c('gtf', 'gtf').", call. = FALSE)
                
                if (!file.exists(annotation_file_query))
                        stop(
                                "The file '",
                                annotation_file_query,
                                " does not seem to exist. Please check the path to your annotation file.",
                                call. = FALSE
                        )
                
                if (!file.exists(annotation_file_subject))
                        stop(
                                "The file '",
                                annotation_file_subject,
                                " does not seem to exist. Please check the path to your annotation file.",
                                call. = FALSE
                        )
                
                message("Start ortholog selection process by ",collapse_by, " ...")
                
                message("Import query annotation file '", annotation_file_query, "' using format ",format[1], " ...")
                
                # import annotation for query
                imported_annotation_qry <-
                        tibble::as_tibble(rtracklayer::import.gff(annotation_file_query, format = format[1]))
                message("Select orthologs in query species ...")
                
                if (format[1] == "gff") {
                        Name <- NULL
                        extracted_features_qry <-
                                dplyr::do(dplyr::group_by(imported_annotation_qry, Name),
                                          extract_features(. , format = format[1]))
                        colnames(extracted_features_qry)[1] <-
                                "query_Name"
                }
                
                if (format[1] == "gtf") {
                        gene_id <- NULL
                        extracted_features_qry <-
                                dplyr::do(dplyr::group_by(imported_annotation_qry, gene_id),
                                          extract_features(. , format = format[1]))
                        colnames(extracted_features_qry)[1] <-
                                "query_Name"
                }
                
                message("Import subject annotation file '", annotation_file_subject, "' using format ",format[2], " ...")
                # import annotation for subjecy
                imported_annotation_sbj <-
                        tibble::as_tibble(rtracklayer::import.gff(annotation_file_subject, format = format[2]))
                
                message("Select orthologs in subject species ...")
                
                if (format[2] == "gff") {
                        Name <- NULL
                        extracted_features_sbj <-
                                dplyr::do(dplyr::group_by(imported_annotation_sbj, Name),
                                          extract_features(. , format = format[2]))
                        colnames(extracted_features_sbj)[1] <-
                                "subject_Name"
                        colnames(extracted_features_sbj)[3] <-
                                "subject_id"
                }
                
                if (format[2] == "gtf") {
                        gene_id <- NULL
                        extracted_features_sbj <-
                                dplyr::do(dplyr::group_by(imported_annotation_sbj, gene_id),
                                          extract_features(. , format = format[2]))
                        colnames(extracted_features_sbj)[1] <-
                                "subject_gene_id"
                        colnames(extracted_features_sbj)[3] <-
                                "subject_id"
                }
                
                message("Join dnds and annotation tables ...")
                # join dnds tbl with query annotation
                joined_dnds_annotation_tbl_qry <-
                        dplyr::inner_join(dplyr::ungroup(dnds_tbl), dplyr::ungroup(extracted_features_qry), by = "query_id")
                
                # join dnds tbl with subject annotation
                joined_dnds_annotation_tbl_sbj <-
                        dplyr::inner_join(dplyr::ungroup(dnds_tbl), dplyr::ungroup(extracted_features_sbj), by = "subject_id")
                # join query and subject tables
                final_join_dnds_annotation_tbl <-
                        dplyr::inner_join(dplyr::ungroup(joined_dnds_annotation_tbl_qry),
                                          dplyr::ungroup(joined_dnds_annotation_tbl_sbj),
                                          by = c("query_id", "subject_id"))
                
                if (collapse_by == "gene_locus") {
                        gene_locus_id <- NULL
                        res <-
                                dplyr::do(
                                        dplyr::group_by(
                                                final_join_dnds_annotation_tbl,
                                                gene_locus_id
                                        ),
                                        filter_best_hits(.)
                                )
                }
                
                if (collapse_by == "splice_variant") {
                        res <- final_join_dnds_annotation_tbl
                }
                
                message("... orthology inference finished successfully!")
                return(res)
        }





