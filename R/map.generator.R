#' @title Generate dNdS Maps Between a Query Organism and Multiple Subject Organisms
#' @description This function allows you to compute dNdS or DS Maps between a query organism
#' and a set subject organisms stored in the same folder. The corresponding dNdS/DS Maps are then stored in an output folder.
#' @param query_file a character string specifying the path to the CDS file of the query organism.
#' @param subjects_folder a character string specifying the path to the folder where CDS files of the subject organisms are stored.
#' @param output_folder a character string specifying the path to the folder where output dNdS/DS Maps should be stored.
#' @param eval a character string specifying the e-value for BLAST based Orthology Inference that is performed
#' in the process of dNdS computations. Please use the scientific notation.
#' @param min_qry_coverage_hsp minimum \code{qcovhsp} (= query coverage of the HSP) of an orthologous hit (a value between 1 and 100).
#' @param min_qry_perc_identity minimum \code{perc_identity} (= percent sequence identity between query and selected HSP) of an orthologous hit (a value between 1 and 100).
#' @param ortho_detection a character string specifying the Orthology Inference method that shall be used to perform
#' dNdS computations. Possible options are: \code{ortho_detection} = \code{"BH"} (BLAST best hit), 
#' \code{ortho_detection} = \code{"RBH"} (BLAST best reciprocal hit), etc.
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type = "multiple"} or \code{aa_aln_type = "pairwise"}. 
#' Default is \code{aa_aln_type = "pairwise"}.
#' @param aa_aln_tool a character string specifying the program that should be used e.g. "clustalw".
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used. 
#' Default is \code{codon_aln_tool = "pal2nal"}. Right now only "pal2nal" can be selected as codon alignment tool.
#' @param dnds_est_method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param comp_cores number of computing cores that shall be used to perform parallelized computations. 
#' @param progress_bar should a progress bar be shown. Default is \code{progress_bar = TRUE}.
#' @param sep a file separator that is used to store maps as csv file.
#' @param ... additional parameters that shall be passed to  \code{\link{dNdS}}.
#' @details
#' Given a query organism and a set of subject organsisms that are stored in the same folder,
#' this function crawls through all subject organsism and as a first step computes the pairwise 
#' dNdS Maps between query and subject organsim and as a second step stores the corresponding Map 
#' in an output folder.
#' @author Hajk-Georg Drost
#' @examples
#' \dontrun{
#' map_generator_dnds(
#'    query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'    subjects_folder = system.file('seqs/map_gen_example', package = 'orthologr'),
#'    aa_aln_type      = "pairwise",
#'    aa_aln_tool      = "NW", 
#'    codon_aln_tool   = "pal2nal", 
#'    dnds_est_method  = "Comeron",
#'    output_folder   = getwd(),
#'    quiet           = TRUE,
#'    comp_cores      = 1
#' )
#' 
#' }
#' @export

map_generator_dnds <- function(query_file, 
                           subjects_folder,
                           output_folder, 
                           eval             = "1E-5",
                           min_qry_coverage_hsp = 50,
                           min_qry_perc_identity = 30,
                           ortho_detection  = "RBH",
                           aa_aln_type      = "pairwise",
                           aa_aln_tool      = "NW", 
                           codon_aln_tool   = "pal2nal", 
                           dnds_est_method  = "Comeron", 
                           comp_cores       = 1,
                           progress_bar     = TRUE,
                           sep              = ";",
                          ... ){
        # retrieve all subject files within a given folder
        subj.files <- list.files(subjects_folder)
        
        if (length(subj.files) == 0)
                stop("Your subject.folder ", subjects_folder, " seems to be empty...", call. = FALSE)
        
        if (!dplyr::between(min_qry_coverage_hsp, 1, 100))
                stop("Please provide a valid min_qry_coverage_hsp value between 1 and 100.", call. = FALSE)
        
        if (!dplyr::between(min_qry_perc_identity, 1, 100))
                stop("Please provide a valid min_qry_perc_identity value between 1 and 100.", call. = FALSE)
        
        message("Starting pairwise genome comparisons (orthology inference and dNdS estimation) between query species: ", basename(query_file), " and subject species: ", paste0(subj.files, collapse = ", "))
        
        # initialize progress bar
        if (progress_bar & (length(subj.files) > 1))
                pb <- utils::txtProgressBar(1, length(subj.files), style = 3)
        
        
        if (!file.exists(output_folder))
                dir.create(output_folder)
        
        qcovhsp <- perc_identity <- NULL
        
        for (i in 1:length(subj.files)) {
                # compute pairwise KaKs/divergence maps between query and all subject files
                OrgQuery_vs_OrgSubj <- dNdS(
                        query_file      = query_file,
                        subject_file    = file.path(subjects_folder, subj.files[i]),
                        eval            = eval,
                        ortho_detection = ortho_detection,
                        aa_aln_type     = aa_aln_type,
                        aa_aln_tool     = aa_aln_tool,
                        dnds_est.method = dnds_est_method,
                        comp_cores      = comp_cores,
                        print_citation = FALSE,
                        ...
                )
                
                message("Filtering for BLAST hits with min_qry_coverage_hsp >= ", min_qry_coverage_hsp, " and min_qry_perc_identity >= ", min_qry_perc_identity, " ...")
                OrgQuery_vs_OrgSubj <- dplyr::filter(OrgQuery_vs_OrgSubj, qcovhsp >= min_qry_coverage_hsp, perc_identity >= min_qry_perc_identity)
                        
                utils::write.table(
                        OrgQuery_vs_OrgSubj,
                        file.path(
                                output_folder,
                                paste0(
                                        "map_q=",
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
        message("\n")
        message("Please cite the following paper when using orthologr for your own research:")
        message("Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. Mol. Biol. Evol. 32 (5): 1221-1231.")
        cat("\n")
        cat(paste0("All maps are stored in ", output_folder, "."))
        
}
        
        


