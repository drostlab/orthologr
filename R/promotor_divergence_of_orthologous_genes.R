#' @title Compute promotor sequence divergence of orthologous genes
#' @description This function computes the promotor sequence divergences of orthologous genes 
#' from a set of pairwise species comparisons. It allows to add this promotor divergence information
#' to an pre-computed \code{dNdS table} generated with \code{\link{dNdS}} or \code{\link{generate_ortholog_tables_all}}.
#' @param promotor_folder a path to a folder containing promotor sequences in \code{fasta} format that
#' were generated with \code{\link[metablastr]{extract_upstream_promotor_seqs}}. 
#' @param ortholog_tables_folder a path to a folder containing dNdS tables generated with \code{\link{generate_ortholog_tables_all}}
#' @param model a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used - must be one of:
#' \itemize{
#' \item \code{K80} (the default)
#' \item \code{raw}
#' \item \code{N}
#' \item \code{TS}
#' \item \code{TV}
#' \item \code{JC69}
#' \item \code{F81} 
#' \item \code{K81}
#' \item \code{F84}
#' \item \code{BH87}
#' \item \code{T92}
#' \item \code{TN93}
#' \item \code{GG95}
#' \item \code{logdet}
#' \item \code{paralin}
#' }
#' @param ortholog_promotor_seq_output a path or name to an output folder where orthologous promotor sequences shall be stored. 
#' @author Hajk-Georg Drost
#' @export

promotor_divergence_of_orthologous_genes <-
        function(promotor_folder,
                 ortholog_tables_folder,
                 model = "K80",
                 ortholog_promotor_seq_output = NULL) {
                
                if (!fs::dir_exists(promotor_folder))
                        stop(
                                "Please provide valid 'promotor_folder' path. The provided directory '",
                                promotor_folder,
                                "' does not seem to exist.",
                                call. = FALSE
                        )
                
                if (!fs::dir_exists(ortholog_tables_folder))
                        stop(
                                "Please provide valid 'ortholog_tables_folder' path. The provided directory '",
                                ortholog_tables_folder,
                                "' does not seem to exist.",
                                call. = FALSE
                        )
                
                if (is.null(ortholog_promotor_seq_output))
                        ortholog_promotor_seq_output <-
                                file.path(getwd(), "ortholog_genes_promotor_seqs")
                
                if (!fs::dir_exists(ortholog_promotor_seq_output)) {
                        message("Creating output folder '",
                                ortholog_promotor_seq_output,
                                "' ...")
                        dir.create(ortholog_promotor_seq_output)
                }
                
                promotor_files <-
                        file.path(promotor_folder, list.files(promotor_folder))
                organism_names <-
                        as.character(sapply(basename(promotor_files), function(x)
                                unlist(stringr::str_split(x, "_"))[1]))
                
                ortholog_tables <-
                        suppressMessages(orthologr::import_ortholog_tables_all(ortholog_tables_folder))
                query_species_name <- unique(ortholog_tables$query_species)
                subject_species_names <- unique(ortholog_tables$subject_species)
                
                organism_names_no_qry <-
                        organism_names[-which(organism_names == query_species_name)]
                promotor_files_no_qry <-
                        promotor_files[-which(organism_names == query_species_name)]
                promotor_file_qry <-
                        promotor_files[which(organism_names == query_species_name)]
                
                message(
                        "Starting promotor sequence extraction of orthologous genes and DNA distance quantification of those promotors for query species: ",
                        query_species_name,
                        " and subject species: ",
                        paste0(subject_species_names, collapse = ", "),
                        " ..."
                )
                
                if (!file.exists(promotor_file_qry))
                        stop(
                                "The promotor file from the query species '",
                                promotor_file_qry,
                                "' does not seem to exist.",
                                call. = FALSE
                        )
                
                promotor_file_qry_sequences <-
                        Biostrings::readDNAStringSet(promotor_file_qry)
                
                promotor_file_qry_sequences@ranges@NAMES <-
                        as.character(sapply(promotor_file_qry_sequences@ranges@NAMES, function(x)
                                unlist(stringr::str_split(x, "_"))[1]))
                
                res <- vector("list", length = length(subject_species_names))
                
                if (!identical(organism_names_no_qry, subject_species_names))
                        stop(
                                "The species names in the 'promotor_folder' and the species names in the 'ortholog_tables_folder' do not match. Please provide
         files for the same species in both folders.",
                                call. = FALSE
                        )
                
                subject_species <- NULL
                
                for (i in seq_len(length(subject_species_names))) {
                        if (organism_names_no_qry[i] != subject_species_names[i])
                                stop(
                                        "The species file used from promotor_folder named '",
                                        organism_names_no_qry[i],
                                        " does not match with the species file used from ortholog_tables_folder named '",
                                        subject_species_names[i],
                                        call. = FALSE
                                )
                        
                        message(
                                "Processing query species: ",
                                query_species_name,
                                " and subject species: ",
                                subject_species_names[i],
                                " ..."
                        )
                        species_specific_ortholog_tables <-
                                dplyr::filter(ortholog_tables,
                                              subject_species == subject_species_names[i])
                        
                        if (!stringr::str_detect(promotor_files_no_qry[i], subject_species_names[i]))
                                stop(
                                        "The 'promotor_file' '",
                                        promotor_files_no_qry[i],
                                        " and the subject_species name from the ortholog_tables '",
                                        subject_species_names[i],
                                        "' do not match.",
                                        call. = FALSE
                                )
                        
                        if (!file.exists(promotor_files_no_qry[i]))
                                stop(
                                        "The promotor file from the subject species '",
                                        promotor_files_no_qry[i],
                                        "' does not seem to exist.",
                                        call. = FALSE
                                )
                        
                        promotor_file_sbj_sequences <-
                                Biostrings::readDNAStringSet(promotor_files_no_qry[i])
                        promotor_file_sbj_sequences@ranges@NAMES <-
                                as.character(sapply(promotor_file_sbj_sequences@ranges@NAMES, function(x)
                                        unlist(stringr::str_split(x, "_"))[1]))
                        
                        qry_matched_genes <-
                                promotor_file_qry_sequences@ranges@NAMES[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables$query_gene_locus_id,
                                                promotor_file_qry_sequences@ranges@NAMES
                                        )
                                )]
                        
                        if (length(qry_matched_genes) == 0)
                                stop(
                                        "The ortholog_table for species '",
                                        subject_species_names[i],
                                        "' does not seem to contain gene_ids that match with the gene_ids in the promotor sequence file '",
                                        promotor_files_no_qry[i],
                                        "'.",
                                        call. = FALSE
                                )
                        
                        message(
                                "Storing promotor sequences of orthologous genes at '",
                                ortholog_promotor_seq_output,
                                "' ..."
                        )
                        qry_output_file_path <-
                                file.path(
                                        ortholog_promotor_seq_output,
                                        paste0(
                                                "ortholog_promotor_seqs_ortholog_",
                                                subject_species_names[i],
                                                "_",
                                                basename(promotor_file_qry)
                                        )
                                )
                        sbj_output_file_path <-
                                file.path(ortholog_promotor_seq_output,
                                          paste0(
                                                  "ortholog_promotor_seqs_",
                                                  basename(promotor_files_no_qry[i])
                                          ))
                        
                        sbj_matched_genes <-
                                promotor_file_sbj_sequences@ranges@NAMES[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables$subject_gene_locus_id,
                                                promotor_file_sbj_sequences@ranges@NAMES
                                        )
                                )]
                        
                        if (length(sbj_matched_genes) == 0)
                                stop(
                                        "The ortholog_table for species '",
                                        query_species_name,
                                        "' does not seem to contain gene_ids that match with the gene_ids in the promotor sequence file '",
                                        promotor_file_sbj_sequences,
                                        "'.",
                                        call. = FALSE
                                )
                        
                        if (length(qry_matched_genes) != length(sbj_matched_genes)) {
                                message(
                                        "The number of orthologous genes that have corresponding promotor sequences in the query (",
                                        length(qry_matched_genes),
                                        " genes) and in the subject (",
                                        length(sbj_matched_genes),
                                        " genes) do not match! Hence, only the intersecting subset will be used."
                                )
                                
                                query_gene_locus_id <- subject_gene_locus_id <- NULL
                                species_specific_ortholog_tables_subset <-
                                        dplyr::filter(
                                                species_specific_ortholog_tables,
                                                query_gene_locus_id %in% qry_matched_genes,
                                                subject_gene_locus_id %in% sbj_matched_genes
                                        )
                                
                                if (nrow(species_specific_ortholog_tables_subset) == 0)
                                        stop(
                                                "Something went wrong. After searching for the subset orhtologous genes, no genes were retained.",
                                                call. = FALSE
                                        )
                                
                                message(
                                        "After subsetting ",
                                        nrow(species_specific_ortholog_tables_subset),
                                        " orthologous genes were retained."
                                )
                                
                                Biostrings::writeXStringSet(promotor_file_qry_sequences[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables_subset$query_gene_locus_id,
                                                promotor_file_qry_sequences@ranges@NAMES
                                        )
                                )], qry_output_file_path)
                                
                                Biostrings::writeXStringSet(promotor_file_sbj_sequences[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables_subset$subject_gene_locus_id,
                                                promotor_file_sbj_sequences@ranges@NAMES
                                        )
                                )], sbj_output_file_path)
                        } else {
                                Biostrings::writeXStringSet(promotor_file_qry_sequences[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables$query_gene_locus_id,
                                                promotor_file_qry_sequences@ranges@NAMES
                                        )
                                )], qry_output_file_path)
                                
                                Biostrings::writeXStringSet(promotor_file_sbj_sequences[stats::na.omit(
                                        match(
                                                species_specific_ortholog_tables$subject_gene_locus_id,
                                                promotor_file_sbj_sequences@ranges@NAMES
                                        )
                                )], sbj_output_file_path)
                        }
                        
                        ortho_promotor_divergence <-
                                promotor_divergence_estimation(query = qry_output_file_path,
                                                                           subject = sbj_output_file_path,
                                                                           model = model)
                        if (nrow(ortho_promotor_divergence) == 0)
                                stop(
                                        "There were no results obtained from the 'promotor_divergence_estimation()' analysis. Please check what might have gone wrong.",
                                        call. = FALSE
                                )
                        
                        res[i] <-
                                list(
                                        dplyr::inner_join(
                                                species_specific_ortholog_tables,
                                                ortho_promotor_divergence,
                                                by = c("query_gene_locus_id", "subject_gene_locus_id")
                                        )
                                )
                        cat("\n")
                }
                
                message("Finished analysis!")
                return(dplyr::bind_rows(res))
                
        }