#' @title A line plot visualizing the number of pairwise orthologs within a ortho table generated with \code{\link{generate_ortholog_tables_all}}
#' @description Given a ortho table generated with \code{\link{generate_ortholog_tables_all}},
#' this function will visualize the number of pairwise orthologs inferred between a reference species A vs a set of subject species B_1, B_2, ...,B_N.
#' @param ortho_tables a \code{ortho tables} that was generated with \code{\link{generate_ortholog_tables_all}}.
#' @param species_order a character string specifying species names listed in the order of phylogenetic/taxonomic distance from the query species.
#' The species names must match with the species names present in the \code{ortho_tables}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param title plot title.
#' @param n_core_orthologs number of core orthologs within the \code{ortho tables}.
#' This number can be retrieved using \code{\link{retrieve_core_orthologs}}.
#' @author Hajk-Georg Drost
#' @export
plot_pairwise_orthologs <- function(ortho_tables, species_order, xlab = "Subject Species", ylab = "Number of reciprocal best hit orthologs", title = "", n_core_orthologs = NULL) {
        
        subject_species <- NULL 
        ortholog_tbl_athaliana_pairwise <- dplyr::summarize(dplyr::group_by(ortho_tables, subject_species),
                                                            n_orthologs = dplyr::n())
        
        ortholog_tbl_athaliana_pairwise$subject_species <- factor(
                ortholog_tbl_athaliana_pairwise$subject_species,
                levels = species_order
        )
        
        subject_species <- n_orthologs <- NULL
        p <- ggplot2::ggplot(
                ortholog_tbl_athaliana_pairwise,
                ggplot2::aes(x = subject_species,
                             y = n_orthologs,
                             group = 1)
        ) + ggplot2::geom_line(size = 2) + ggplot2::geom_point(size = 4) + ggplot2::geom_abline(
                intercept = ifelse(!is.null(n_core_orthologs), n_core_orthologs, 0),
                size = 2,
                col = "darkred",
                alpha = 0.4
        ) +
                ggplot2::geom_text(
                        ggplot2::aes(label = n_orthologs),
                        hjust = 0,
                        vjust = -1.5,
                        size = 6
                ) + ggplot2::scale_y_continuous(limits = c(0, max(ortholog_tbl_athaliana_pairwise$n_orthologs) + (max(ortholog_tbl_athaliana_pairwise$n_orthologs) / 10)),
                                                breaks = scales::pretty_breaks(n = 10)) +
                ggplot2::theme_minimal() +
                ggplot2::labs(x = xlab, y = ylab, title = title) +
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
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
        
        
        return(p)
}

