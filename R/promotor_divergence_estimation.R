#' @title Estimate the DNA distance between promotor sequences
#' @description Given two \code{fasta} files storing the promotor sequence
#' of protein coding genes, this function estimates the pairwise DNA distance between these promotor sequences.
#' @param query file path to a query \code{fasta} file storing promotor sequences of interest (e.g. generated with \code{\link[metablastr]{extract_upstream_promotor_seqs}}).
#' @param subject file path to a subject \code{fasta} file storing promotor sequences of interest (e.g. generated with \code{\link[metablastr]{extract_upstream_promotor_seqs}}).
#' @param model a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used - must be one of:
#' \itemize{
#' \item  \code{K80} (the default)
#' \item \code{raw}
#' \item  \code{N}
#' \item  \code{TS}
#' \item  \code{TV}
#' \item  \code{JC69}
#' \item  \code{F81} 
#' \item \code{K81}
#' \item \code{F84}
#' \item \code{BH87}
#' \item \code{T92}
#' \item \code{TN93}
#' \item \code{GG95}
#' \item \code{logdet}
#' \item \code{paralin}
#' }
#' @author Hajk-Georg Drost
#' @export

promotor_divergence_estimation <-
        function(query,
                 subject,
                 model = "K80") {
                
                if (!file.exists(query))
                        stop("Please provide a valid path to a query file.", call. = FALSE)
                
                if (!file.exists(subject))
                        stop("Please provide a valid path to a subject file.", call. = FALSE)
                
                message("Starting promotor divergence estimation between '", query, "' and '", subject, "' ...")
                
                # import promotor sequences
                query_seqs <- Biostrings::readDNAStringSet(query)
                subject_seqs <- Biostrings::readDNAStringSet(subject)
                
                if (length(query_seqs) != length(subject_seqs))
                        stop("The files '",query,"' and '",subject,"' don't contain the same number of promotor sequences. Please provide query and subject files that contain the same number of promotor sequences.", call. = FALSE)
                
                tmp_file_name_qry <- unlist(stringr::str_split(basename(query), "[.]"))[1]
                tmp_file_name_sbj <- unlist(stringr::str_split(basename(subject), "[.]"))[1]
                
                # compute global pairwise alignments between query and subject promotor sequences
                Alignments <-
                        Biostrings::pairwiseAlignment(query_seqs, subject_seqs)
                
                names_query_seqs <- names(query_seqs)
                names_subject_seqs <- names(subject_seqs)
                
                # retrieve query alignment sequences
                qry_seqs_with_gaps <-
                        Biostrings::BStringSet(Alignments@pattern)
                names(qry_seqs_with_gaps) <- names(names_query_seqs)
                
                # save query part of alignment sequences for import into the ape package
                Biostrings::writeXStringSet(
                        qry_seqs_with_gaps,
                        file.path(tempdir(), paste0(tmp_file_name_qry, "_query_seqs_with_gaps.fasta")) 
                )
                
                # import query part of alignment into ape
                qry_seqs_with_gaps_ape <-
                        ape::read.dna(file.path(tempdir(), paste0(tmp_file_name_qry, "_query_seqs_with_gaps.fasta")),
                                      format = "fasta")
                
                # retrieve subject part of alignment sequences
                sbj_seqs_with_gaps <-
                        Biostrings::BStringSet(Alignments@subject)
                names(sbj_seqs_with_gaps) <- names(names_subject_seqs)
                
                # save subject part of alignment sequences for import into the ape package
                Biostrings::writeXStringSet(
                        sbj_seqs_with_gaps,
                        file.path(tempdir(), paste0(tmp_file_name_sbj, "_subject_seqs_with_gaps.fasta"))
                )
                
                # import subject part of alignment sequences into ape
                sbj_seqs_with_gaps_ape <-
                        ape::read.dna(file.path(tempdir(), paste0(tmp_file_name_sbj, "_subject_seqs_with_gaps.fasta")),
                                      format = "fasta")
                
                dist_vals <- vector("numeric", length(sbj_seqs_with_gaps_ape))
                indel_blocks <- vector("numeric", length(sbj_seqs_with_gaps_ape))
                indel <- vector("numeric", length(sbj_seqs_with_gaps_ape))
                
                for (i in seq_len(length(sbj_seqs_with_gaps_ape))) {
                        
                        ape::write.FASTA(
                                sbj_seqs_with_gaps_ape[i],
                                file.path(tempdir(), paste0(tmp_file_name_sbj, "_seq_",i, ".fasta"))
                        )
                        
                        ape::write.FASTA(
                                qry_seqs_with_gaps_ape[i],
                                file.path(tempdir(), paste0(tmp_file_name_sbj, "_seq_",i, ".fasta")),
                                append = TRUE
                        )
                        
                        seq_i <-
                                ape::read.dna(file.path(tempdir(), paste0(tmp_file_name_sbj, "_seq_",i, ".fasta")), format = "fasta")
                        dist_vals[i] <- ape::dist.dna(seq_i, model = model, variance = TRUE)
                        indel_blocks[i] <- ape::dist.dna(seq_i, model = "indelblock")
                        indel[i] <- ape::dist.dna(seq_i, model = "indel")
                }
                
                res <- tibble::tibble(
                        query_gene_locus_id = names(query_seqs),
                        subject_gene_locus_id = names(subject_seqs),
                        promotor_name = paste0(names(query_seqs), "_", names(subject_seqs)),
                        promotor_dist = dist_vals,
                        promotor_aln_score = Alignments@score,
                        promotor_indel_blocks = indel_blocks,
                        promotor_indels = indel
                )
                
                return(res)
        }
