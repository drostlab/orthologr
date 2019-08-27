#' @title Diverse line plots visualizing the number of pairwise orthologs within a ortho table generated with \code{\link{generate_ortholog_tables_all}}
#' based on different sets of homology thresholds.
#' @description Given a ortho table generated with \code{\link{generate_ortholog_tables_all}},
#' this function will visualize the number of pairwise orthologs inferred between a reference species A vs a set of subject species B_1, B_2, ...,B_N
#' using a diverse range of homology thresholds to enable an analytical decision for chosing homology thresholds to define orthologous genes.
#' @param ortho_tables a \code{ortho tables} that was generated with \code{\link{generate_ortholog_tables_all}}.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in the \code{ortho_tables}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title plot title.
#' @author Hajk-Georg Drost
#' @export
plot_diverse_homology_thresholds <-
        function(ortho_tables,
                 species_order,
                 xlab = "Subject Species",
                 ylab = "Number of reciprocal best hit orthologs",
                 title = "") {
                qcovhsp <-
                        scope <- perc_identity <- subject_species <- n_orthologs <- NULL
                
                qcovhsp_70_perc_identity_30 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 70,
                                        scope >= 0.7,
                                        perc_identity >= 30
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_70_perc_identity_50 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 70,
                                        scope >= 0.7,
                                        perc_identity >= 50
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_70_perc_identity_70 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 70,
                                        scope >= 0.7,
                                        perc_identity >= 70
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_70_perc_identity_90 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 70,
                                        scope >= 0.7,
                                        perc_identity >= 90
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_30_perc_identity_30 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 30,
                                        scope >= 0.3,
                                        perc_identity >= 30
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_50_perc_identity_30 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 50,
                                        scope >= 0.5,
                                        perc_identity >= 30
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_90_perc_identity_30 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 90,
                                        scope >= 0.9,
                                        perc_identity >= 30
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_90_perc_identity_50 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 90,
                                        scope >= 0.9,
                                        perc_identity >= 50
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_90_perc_identity_70 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 90,
                                        scope >= 0.9,
                                        perc_identity >= 70
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                qcovhsp_90_perc_identity_90 <-
                        dplyr::summarize(dplyr::group_by(
                                dplyr::filter(
                                        ortho_tables,
                                        qcovhsp >= 90,
                                        scope >= 0.9,
                                        perc_identity >= 90
                                ),
                                subject_species
                        ),
                        n_orthologs = dplyr::n())
                
                p_all_ortho_thresholds <-
                        orthologr::plot_pairwise_orthologs(
                                ortho_tables = ortho_tables,
                                species_order = species_order,
                                n_core_orthologs = NULL,
                                xlab = xlab,
                                ylab = ylab,
                                title = title
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_70_perc_identity_30,
                                size = 2,
                                colour = "#3B4992FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_70_perc_identity_30,
                                                 colour = "#3B4992FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_70_perc_identity_30,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#3B4992FF"
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_70_perc_identity_50,
                                size = 2,
                                colour = "#808180FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_70_perc_identity_50,
                                                 colour = "#808180FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_70_perc_identity_50,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#808180FF"
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_70_perc_identity_70,
                                size = 2,
                                colour = "#EE0000FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_70_perc_identity_70,
                                                 colour = "#EE0000FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_70_perc_identity_70,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#EE0000FF"
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_70_perc_identity_90,
                                size = 2,
                                colour = "#008B45FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_70_perc_identity_90,
                                                 colour = "#008B45FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_70_perc_identity_90,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#008B45FF"
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_30_perc_identity_30,
                                size = 2,
                                colour = "#631879FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_30_perc_identity_30,
                                                 colour = "#631879FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_30_perc_identity_30,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#631879FF"
                        ) + ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_50_perc_identity_30,
                                size = 2,
                                colour = "#008280FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_50_perc_identity_30,
                                                 colour = "#008280FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_50_perc_identity_30,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#008280FF"
                        ) +
                        ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_90_perc_identity_30,
                                size = 2,
                                colour = "#BB0021FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_90_perc_identity_30,
                                                 colour = "#BB0021FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_90_perc_identity_30,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#BB0021FF"
                        ) +
                        ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_90_perc_identity_50,
                                size = 2,
                                colour = "#5F559BFF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_90_perc_identity_50,
                                                 colour = "#5F559BFF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_90_perc_identity_50,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#5F559BFF"
                        ) +
                        ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_90_perc_identity_70,
                                size = 2,
                                colour = "darkgreen"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_90_perc_identity_70,
                                                 colour = "darkgreen") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_90_perc_identity_70,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "darkgreen"
                        ) +
                        ggplot2::geom_line(
                                ggplot2::aes(
                                        x = subject_species,
                                        y = n_orthologs,
                                        group = 1
                                ),
                                data = qcovhsp_90_perc_identity_90,
                                size = 2,
                                colour = "#A20056FF"
                        )  + ggplot2::geom_point(size = 4,
                                                 data = qcovhsp_90_perc_identity_90,
                                                 colour = "#A20056FF") +
                        ggplot2::geom_text(
                                ggplot2::aes(label = n_orthologs),
                                data = qcovhsp_90_perc_identity_90,
                                hjust = 0,
                                vjust = 3,
                                size = 4,
                                colour = "#A20056FF"
                        ) + ggplot2::scale_fill_manual(
                                values = c(
                                        "#3B4992FF",
                                        "#808180FF",
                                        "#EE0000FF",
                                        "#008B45FF",
                                        "#631879FF",
                                        "#008280FF",
                                        "#BB0021FF",
                                        "#5F559BFF",
                                        "darkgreen",
                                        "#A20056FF"
                                ),
                                name = "Parameter Options",
                                labels = c(
                                        "qcovhsp_70_perc_identity_30",
                                        "qcovhsp_70_perc_identity_50",
                                        "qcovhsp_70_perc_identity_70",
                                        "qcovhsp_70_perc_identity_90",
                                        "qcovhsp_30_perc_identity_30",
                                        "qcovhsp_50_perc_identity_30",
                                        "qcovhsp_90_perc_identity_30",
                                        "qcovhsp_90_perc_identity_50",
                                        "qcovhsp_90_perc_identity_70",
                                        "qcovhsp_90_perc_identity_90"
                                )
                        )
                
                return(p_all_ortho_thresholds)
        }