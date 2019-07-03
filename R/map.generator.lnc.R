#' @title Infer orthologous lncRNAs between muliple species
#' @description Inference of orthologous lncRNAs between multiple species is performed via pairwise comparisons.
#'  The corresponding orthologous tables are then stored in an output folder.
#' @param query_file a character string specifying the path to the lncRNAs file of the query organism in \code{fasta} format.
#' @param subjects_folder a character string specifying the path to the folder where lncRNAs files in \code{fasta} format of the subject organisms are stored.
#' @param output_folder a character string specifying the path to the folder where output orthologous tables should be stored.
#' @param eval a character string specifying the e-value for BLAST based orthology inference. Please use the scientific notation.
#' @param ortho_detection a character string specifying the Orthology Inference method that shall be used to perform
#' dNdS computations. Possible options are: 
#' \itemize{
#'  \item \code{ortho_detection = "BH"}: BLAST best unidirectional hit
#'  \item \code{ortho_detection = "RBH"}: BLAST best reciprocal hit
#' }
#' @param min_qry_coverage_hsp minimum \code{qcovhsp} (= query coverage of the HSP) of an orthologous hit (a value between 1 and 100).
#' @param min_qry_perc_identity minimum \code{perc_identity} (= percent sequence identity between query and selected HSP) of an orthologous hit (a value between 1 and 100).
#' @param comp_cores number of computing cores that shall be used to perform parallelized computations. 
#' @param progress_bar should a progress bar be shown. Default is \code{progress_bar = TRUE}.
#' @param sep a file separator that is used to store maps as csv file.
#' @param ... additional parameters that shall be passed to  \code{\link{dNdS}}.
#' @details
#' Given a query organism and a set of subject organsisms that are stored in the same folder,
#' this function crawls through all subject organsism and infers the lncRNA homologs in
#' pairwise species comparisons.
#' @author Hajk-Georg Drost
#' @examples
#' \dontrun{
#' map_generator_lnc(
#'    query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'    subjects_folder = system.file('seqs/map_gen_example', package = 'orthologr'),
#'    output_folder   = getwd(),
#'    quiet           = TRUE,
#'    comp_cores      = 1
#' )
#' }
#' @export

map_generator_lnc <- function(query_file, 
                          subjects_folder,
                          output_folder, 
                          eval             = "1E-5",
                          ortho_detection  = "RBH",
                          min_qry_coverage_hsp = 50,
                          min_qry_perc_identity = 30,
                          comp_cores       = 1,
                          progress_bar     = TRUE,
                          sep              = ";",
                          ... ){
        
        if (!fs::file_exists(query_file))
                stop("Please provide a valid path to the 'query_file'.", call. = FALSE)
        
        # retrieve all subject files within a given folder
        subj.files <- list.files(subjects_folder)
        
        if (length(subj.files) == 0)
                stop("Your subject.folder ", subjects_folder, " seems to be empty...", call. = FALSE)
        
        message("Starting pairwise orthology inference of lncRNAs between query species: ", basename(query_file), " and subject species: ", paste0(subj.files, collapse = ", "))
        
        # initialize progress bar
        if (progress_bar & (length(subj.files) > 1))
                pb <- utils::txtProgressBar(1, length(subj.files), style = 3)
        
        if (!file.exists(output_folder))
                dir.create(output_folder)
        
        qcovhsp <- perc_identity <- NULL
        
        for (i in seq_len(length(subj.files))) {
                
                message("LncRNA orthology inference between ", basename(query_file), " and ", subj.files[i], " (",i ,"/", length(subj.files),")")
                # perform pairwise lncRNA orthology inference
                OrgQuery_vs_OrgSubj <- orthologs_lnc(
                        query_file      = query_file,
                        subject_file    = file.path(subjects_folder, subj.files[i]),
                        eval            = eval,
                        ortho_detection = ortho_detection,
                        comp_cores      = comp_cores,
                        ...
                )
                
                message("Filtering for BLAST hits with min_qry_coverage_hsp >= ", min_qry_coverage_hsp, " and min_qry_perc_identity >= ", min_qry_perc_identity, " ...")
                OrgQuery_vs_OrgSubj <- dplyr::filter(OrgQuery_vs_OrgSubj, qcovhsp >= min_qry_coverage_hsp, perc_identity >= min_qry_perc_identity)
                
                utils::write.table(
                        OrgQuery_vs_OrgSubj,
                        file.path(
                                output_folder,
                                paste0(
                                        "lncRNA_orthologs_q=",
                                        basename(query_file),
                                        "_s=",
                                        subj.files[i],
                                        "_orthodetec=",
                                        ortho_detection,
                                        "_eval=",
                                        eval,
                                        ".csv"
                                )
                        ),
                        sep       = sep,
                        col.names = TRUE,
                        row.names = FALSE,
                        quote     = FALSE
                )
                
                if (progress_bar)
                        utils::setTxtProgressBar(pb, i)
                
        }
        message("\n")
        message("All maps are stored in ", output_folder, ".")
        message("Orthology inference finished successfully!")
}




