#' @title Select orthologs based on either gene locus or splice variant
#' @description This function selects orthologs based on either gene locus or splice variant.
#' @param dnds_tbl a tibble returned by \code{link{dNdS}}.
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
                 format = c("gtf", "gtf")) {
                if (!is.element(collapse_by, c("gene_locus", "splice_variant")))
                        stop(
                                "Please specify a collapse method provided by this function: collapse_by = 'gene_locus' or collapse_by = 'splice_variant'.",
                                call. = FALSE
                        )
                
                if (!all(is.element(format, c("gtf", "gff"))))
                        stop(
                                "Please provide a vector of length 2 that contains available formats: either 'gtf' or 'gff'.",
                                call. = FALSE
                        )
                
                if (length(format) != 2)
                        stop(
                                "Please provide a vector of length 2 to specify the query and subject annotation file formats. E.g. format = c('gtf', 'gtf').",
                                call. = FALSE
                        )
                
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
                
                temp_store <-
                        file.path(tempdir(),
                                  paste0(
                                          "final_join_dnds_annotation_tbl_",
                                          basename(annotation_file_query)
                                  ))
                
                if (file.exists(temp_store)) {
                        message("Import pre-computed file final_join_dnds_annotation_tbl ...")
                        final_join_dnds_annotation_tbl <-
                                readRDS(temp_store)
                        
                } else {
                        message("Start ortholog selection process by ",
                                collapse_by,
                                " ...")
                        
                        message(
                                "Import query annotation file '",
                                annotation_file_query,
                                "' using format ",
                                format[1],
                                " ..."
                        )
                        
                        # import annotation for query
                        imported_annotation_qry <-
                                tibble::as_tibble(rtracklayer::import.gff(annotation_file_query, format = format[1]))
                        
                        message("Select orthologs in query species ...")
                        
                        if (format[1] == "gff") {
                                Name <- NULL
                                extracted_features_qry <-
                                        dplyr::do(
                                                dplyr::group_by(imported_annotation_qry, Name),
                                                extract_features(. , format = format[1])
                                        )
                                colnames(extracted_features_qry)[2] <-
                                        "query_gene_locus_id"
                                
                                extracted_features_qry <-
                                        dplyr::select(
                                                dplyr::ungroup(extracted_features_qry),
                                                "query_gene_locus_id",
                                                "query_id"
                                        )
                                
                        }
                        
                        if (format[1] == "gtf") {
                                gene_id <- NULL
                                extracted_features_qry <-
                                        dplyr::do(
                                                dplyr::group_by(imported_annotation_qry, gene_id),
                                                extract_features(. , format = format[1])
                                        )
                                colnames(extracted_features_qry)[2] <-
                                        "query_gene_locus_id"
                                
                                extracted_features_qry <-
                                        dplyr::select(
                                                dplyr::ungroup(extracted_features_qry),
                                                "query_gene_locus_id",
                                                "query_id"
                                        )
                        }
                        
                        message(
                                "Import subject annotation file '",
                                annotation_file_subject,
                                "' using format ",
                                format[2],
                                " ..."
                        )
                        # import annotation for subjecy
                        imported_annotation_sbj <-
                                tibble::as_tibble(
                                        rtracklayer::import.gff(annotation_file_subject, format = format[2])
                                )
                        
                        message("Select orthologs in subject species ...")
                        
                        if (format[2] == "gff") {
                                Name <- NULL
                                extracted_features_sbj <-
                                        dplyr::do(
                                                dplyr::group_by(imported_annotation_sbj, Name),
                                                extract_features(. , format = format[2])
                                        )
                                colnames(extracted_features_sbj)[2] <-
                                        "subject_gene_locus_id"
                                colnames(extracted_features_sbj)[3] <-
                                        "subject_id"
                                
                                extracted_features_sbj <-
                                        dplyr::select(
                                                dplyr::ungroup(extracted_features_sbj),
                                                "subject_gene_locus_id",
                                                "subject_id"
                                        )
                                
                        }
                        
                        if (format[2] == "gtf") {
                                gene_id <- NULL
                                extracted_features_sbj <-
                                        dplyr::do(
                                                dplyr::group_by(imported_annotation_sbj, gene_id),
                                                extract_features(. , format = format[2])
                                        )
                                colnames(extracted_features_sbj)[2] <-
                                        "subject_gene_locus_id"
                                colnames(extracted_features_sbj)[3] <-
                                        "subject_id"
                                
                                extracted_features_sbj <-
                                        dplyr::select(
                                                dplyr::ungroup(extracted_features_sbj),
                                                "subject_gene_locus_id",
                                                "subject_id"
                                        )
                        }
                        
                        message("Join dnds and annotation tables ...")
                        # join dnds tbl with query annotation
                        joined_dnds_annotation_tbl_qry <-
                                dplyr::inner_join(dnds_tbl,
                                                  extracted_features_qry,
                                                  by = "query_id")
                        
                        
                        # join dnds tbl with subject annotation
                        joined_dnds_annotation_tbl_sbj <-
                                dplyr::inner_join(dnds_tbl,
                                                  extracted_features_sbj,
                                                  by = "subject_id")
                        
                        joined_dnds_annotation_tbl_sbj <-
                                dplyr::select(
                                        joined_dnds_annotation_tbl_sbj,
                                        "query_id",
                                        "subject_gene_locus_id",
                                        "subject_id"
                                )
                        
                        # join query and subject tables
                        final_join_dnds_annotation_tbl <-
                                dplyr::inner_join(
                                        joined_dnds_annotation_tbl_qry,
                                        joined_dnds_annotation_tbl_sbj,
                                        by = c("query_id", "subject_id")
                                )
                        
                        if (!file.exists(temp_store)) {
                                saveRDS(final_join_dnds_annotation_tbl,
                                        temp_store)
                        }
                }
                
                if (collapse_by == "gene_locus") {
                        message(
                                "Select splice variant with smallest e-value for each gene locus of query species ..."
                        )
                        query_gene_locus_id <-
                                subject_gene_locus_id <- NULL
                        res_qry <-
                                dplyr::do(
                                        dplyr::group_by(
                                                final_join_dnds_annotation_tbl,
                                                query_gene_locus_id
                                        ),
                                        filter_best_hits(.)
                                )
                        
                        message(
                                "Select splice variant with smallest e-value for each gene locus of subject species ..."
                        )
                        res_sbj <-
                                dplyr::ungroup(dplyr::do(
                                        dplyr::group_by(
                                                final_join_dnds_annotation_tbl,
                                                subject_gene_locus_id
                                        ),
                                        filter_best_hits(.)
                                ))
                        
                        res_sbj <-
                                dplyr::select(dplyr::ungroup(res_sbj),
                                              "query_id",
                                              "subject_id")
                        
                        res <-
                                dplyr::inner_join(dplyr::ungroup(res_qry),
                                                  res_sbj,
                                                  by = c("query_id", "subject_id"))
                        
                        if (length(unique(res$query_id)) == length(unique(res$query_gene_locus_id))) {
                                message("Good News: ")
                                message(
                                        "Number of unique query_ids matches the number of unique query_gene_locus_ids. Thus, each gene locus has now a representative splice variant. The one having the smallest evalue and longest alignment length if two splice variants have the same evalue."
                                )
                                message("\n")
                        } else {
                                message(
                                        "Please check your input data. It seems that the number of unique query_ids DOES NOT match the number of unique query_gene_locus_ids. Thus, some query gene loci seem to have multiple splice variants as representatives."
                                )
                                message("\n")
                        }
                        
                        if (length(unique(res$subject_id)) == length(unique(res$subject_gene_locus_id))) {
                                message("Good news:")
                                message(
                                        "Number of unique subject_ids matches the number of unique subject_gene_locus_ids Thus, each gene locus has now a representative splice variant. The one having the smallest evalue and longest alignment length if two splice variants have the same evalue."
                                )
                                message("\n")
                        } else {
                                message(
                                        "Please check your input data. It seems that the number of unique subject_ids DOES NOT match the number of unique subject_gene_locus_ids Thus, some subject gene loci seem to have multiple splice variants as representatives."
                                )
                                message("\n")
                        }
                        
                        if (length(unique(dnds_tbl$query_id)) < length(unique(res$query_id)))
                                stop(
                                        "Something must have gone wrong... The unique number of query_ids in the dNdS_tbl is smaller than the unique number of query_ids in the collapsed output file. Please check your input data.",
                                        call. = FALSE
                                )
                        
                }
                
                if (collapse_by == "splice_variant") {
                        res <- final_join_dnds_annotation_tbl
                        query_id <- subject_id <- NULL
                        res <- dplyr::distinct(
                                res,
                                query_id,
                                query_gene_locus_id,
                                subject_id,
                                subject_gene_locus_id,
                                .keep_all = TRUE
                        )
                }
                
                qry_species_name <-
                        unlist(stringr::str_split(basename(annotation_file_query), "_"))[1]
                sbj_species_name <-
                        unlist(stringr::str_split(basename(annotation_file_subject), "_"))[1]
                res <-
                        dplyr::mutate(
                                res,
                                qry_species = rep(qry_species_name, nrow(res)),
                                subject_species = rep(sbj_species_name, nrow(res))
                        )
                
                evalue <- bit_score <- perc_identity <- NULL
                s_end <-
                        dN <-
                        dS <- qry_species <- subject_species <- NULL
                res <- dplyr::select(
                        res,
                        qry_species,
                        subject_species,
                        query_id,
                        query_gene_locus_id,
                        subject_id,
                        subject_gene_locus_id,
                        dN,
                        dS,
                        dNdS,
                        evalue,
                        bit_score,
                        perc_identity:s_end
                )
                
                message("... orthology inference finished successfully!")
                return(res)
        }
