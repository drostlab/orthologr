#' @title Diverse line plots visualizing the number of core orthologs within a ortho table generated with \code{\link{generate_ortholog_tables_all}}
#' based on different sets of homology thresholds.
#' @description Given a ortho table generated with \code{\link{generate_ortholog_tables_all}},
#' this function will visualize the number of core orthologs
#' using a diverse range of homology thresholds to enable an analytical decision for chosing homology thresholds to define orthologous genes.
#' @param ortho_tables a \code{ortho tables} that was generated with \code{\link{generate_ortholog_tables_all}}.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in the \code{ortho_tables}.
#' @param type type of core orthologs that shall be visualized:
#' \itemize{
#' \item \code{type = "gene_locus"}: visualizes the number of core orthologs for a gene locus (hence, different splice variants can represent the same orthologous gene locus relationship)
#' \item \code{type = "both"}: visualizes the difference in numbers of core orthologs for a gene locus versus orthologous splice variant (thus, only the same splice variant can represent the same orthologous gene locus relationship)
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title plot title.
#' @author Hajk-Georg Drost
#' @export
plot_diverse_homology_thresholds_core_orthologs <-
        function(ortho_tables,
                 species_order,
                 type = "both",
                 xlab = "Parameter set for selecting core orthologs",
                 ylab = "Number of corresponding core orthologs",
                 title = "") {
                if (!is.element(type, c("gene_locus", "both")))
                        stop(
                                "Please provide a valid 'type' argument: type = 'gene_locus' or type = 'both'.",
                                call. = FALSE
                        )
                
                qcovhsp <- perc_identity <- scope <- NULL
                core_sets_multi <-
                        dplyr::bind_rows(
                                list(
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 50,
                                                        perc_identity >= 30,
                                                        scope >= 0.5
                                                ),
                                                "qcovhsp_50_perc_identity_30",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 30,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_30",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 50,
                                                        perc_identity >= 50,
                                                        scope >= 0.5
                                                ),
                                                "qcovhsp_50_perc_identity_50",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 50,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_50",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 70,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_70",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 90,
                                                        perc_identity >= 70,
                                                        scope >= 0.9
                                                ),
                                                "qcovhsp_90_perc_identity_70",
                                                species_order = species_order
                                        ),
                                        testCoreOrthoParamsGeneLocus(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 90,
                                                        perc_identity >= 90,
                                                        scope >= 0.9
                                                ),
                                                "qcovhsp_90_perc_identity_90",
                                                species_order = species_order
                                        )
                                )
                        )
                
                core_sets_multi$parameter <- factor(
                        core_sets_multi$parameter,
                        levels = c(
                                "qcovhsp_50_perc_identity_30",
                                "qcovhsp_70_perc_identity_30",
                                "qcovhsp_50_perc_identity_50",
                                "qcovhsp_70_perc_identity_50",
                                "qcovhsp_70_perc_identity_70",
                                "qcovhsp_90_perc_identity_70",
                                "qcovhsp_90_perc_identity_90"
                        )
                )
                
                parameter <- n_core_genes <- NULL
                p_core_sets_multi <- ggplot2::ggplot(core_sets_multi,
                                                     ggplot2::aes(
                                                             x = parameter,
                                                             y = n_core_genes,
                                                             group = 1
                                                     )) + ggplot2::geom_line(size = 2) + ggplot2::geom_point(size = 4) + ggplot2::geom_text(
                                                             ggplot2::aes(label = n_core_genes),
                                                             hjust = 0,
                                                             vjust = -1.5,
                                                             size = 6
                                                     ) + ggplot2::scale_y_continuous(limits = c(0, max(core_sets_multi$n_core_genes) + (
                                                             max(core_sets_multi$n_core_genes) / 10
                                                     )),
                                                     breaks = scales::pretty_breaks(n = 10)) +
                        ggplot2::theme_minimal() +
                        ggplot2::labs(x = xlab,
                                      y = ylab,
                                      title = title) +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 18, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 18, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 18,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(
                                angle = 90,
                                vjust = 1,
                                hjust = 1
                        ))
                
                core_sets_multsplice_variant <-
                        dplyr::bind_rows(
                                list(
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 50,
                                                        perc_identity >= 30,
                                                        scope >= 0.5
                                                ),
                                                "qcovhsp_50_perc_identity_30",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 30,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_30",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 50,
                                                        perc_identity >= 50,
                                                        scope >= 0.5
                                                ),
                                                "qcovhsp_50_perc_identity_50",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 50,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_50",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 70,
                                                        perc_identity >= 70,
                                                        scope >= 0.7
                                                ),
                                                "qcovhsp_70_perc_identity_70",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 90,
                                                        perc_identity >= 70,
                                                        scope >= 0.9
                                                ),
                                                "qcovhsp_90_perc_identity_70",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        ),
                                        testCoreOrthoParamsSpliceVariant(
                                                dplyr::filter(
                                                        ortho_tables,
                                                        qcovhsp >= 90,
                                                        perc_identity >= 90,
                                                        scope >= 0.9
                                                ),
                                                "qcovhsp_90_perc_identity_90",
                                                species_order = species_order,
                                                n_core_species = length(species_order)
                                        )
                                )
                        )
                
                
                core_sets_multsplice_variant$parameter <- factor(
                        core_sets_multsplice_variant$parameter,
                        levels = c(
                                "qcovhsp_50_perc_identity_30",
                                "qcovhsp_70_perc_identity_30",
                                "qcovhsp_50_perc_identity_50",
                                "qcovhsp_70_perc_identity_50",
                                "qcovhsp_70_perc_identity_70",
                                "qcovhsp_90_perc_identity_70",
                                "qcovhsp_90_perc_identity_90"
                        )
                )
                
                p_core_sets_multsplice_variant <-
                        p_core_sets_multi + ggplot2::geom_line(size = 2,
                                                               data = core_sets_multsplice_variant,
                                                               colour = "#008B45FF") + ggplot2::geom_point(size = 4,
                                                                                                           data = core_sets_multsplice_variant,
                                                                                                           colour = "#008B45FF") + ggplot2::geom_text(
                                                                                                                   ggplot2::aes(label = n_core_genes),
                                                                                                                   data = core_sets_multsplice_variant,
                                                                                                                   colour = "#008B45FF",
                                                                                                                   hjust = 0,
                                                                                                                   vjust = 3,
                                                                                                                   size = 6
                                                                                                           )
                
                if (type == "gene_locus")
                        return(p_core_sets_multi)
                
                if (type == "both")
                        return(p_core_sets_multsplice_variant)
                
        }
