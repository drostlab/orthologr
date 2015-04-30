#' @title Generate dNdS Maps Between a Query Organism and Multiple Subject Organisms
#' @description This function allows you to compute dNdS or DS Maps between a query organism
#' and a set subject organisms stored in the same folder. The corresponding dNdS/DS Maps are then stored in an output folder.
#' @param query_file a character string specifying the path to the CDS file of the query organism.
#' @param subjects.folder a character string specifying the path to the folder where CDS files of the subject organisms are stored.
#' @param output.folder a character string specifying the path to the folder where output dNdS/DS Maps should be stored stored.
#' @param eval a character string specifying the e-value for BLAST based Orthology Inference that is performed
#' in the process of dNdS computations. Please use the scientific notation.
#' @param ortho_detection a character string specifying the Orthology Inference method that shall be used to perform
#' dNdS computations. Possible options are: \code{ortho_detection} = \code{"BH"} (BLAST best hit), 
#' \code{ortho_detection} = \code{"RBH"} (BLAST best reciprocal hit), etc.
#' @param comp_cores number of computing cores that shall be used to perform parallelized computations. 
#' @param ds.values a logical value indicating whether divegrence stratum values (ds values) or dNdS values should be returned by \code{\link{divergence_stratigraphy}}. Default is \code{ds.values = FALSE}.
#' @param progress.bar should a progress bar be shown. Default is \code{progress.bar = TRUE}.
#' @param sep a file separator that is used to store maps as csv file.
#' @param ... additional parameters that shall be passed to  \code{\link{divergence_stratigraphy}}.
#' @details
#' 
#' Given a query organism and a set of subject organsisms that are stored in the same folder,
#' this function crawls through all subject organsism and as a first step computes the pairwise 
#' dNdS or DS Maps between query and subject organsim and as a second step stores the corresponding Map 
#' in an output folder.
#' @author Hajk-Georg Drost
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' map.generator(
#'    query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'    subjects.folder = system.file('seqs/map_gen_example', package = 'orthologr'),
#'    output.folder = getwd(),
#'    quiet = TRUE
#' )
#' 
#' }
#' @export

map.generator <- function(query_file, 
                           subjects.folder,
                           output.folder, 
                           eval             = "1E-5",
                           ortho_detection  = "RBH", 
                           comp_cores       = 1,
                           ds.values         = FALSE,
                           subject.id       = TRUE,
                           progress.bar     = TRUE,
                           sep              = ";",
                          ... ){
        
        # retrieve all subject files within a given folder
        subj.files <- list.files(subjects.folder)
        
        if(length(subj.files) == 0)
                stop("your subject.folder ",subjects.folder, " is empty...")
        
        # initialize progress bar
        if(progress.bar)
                pb   <- txtProgressBar(1, length(subj.files), style = 3)
        
        for(i in 1:length(subj.files)){
                # compute pairwise KaKs/divergence maps between query and all subject files 
                OrgQuery_vs_OrgSubj <- divergence_stratigraphy(
                                              query_file      = query_file,
                                              subject_file    = file.path(subjects.folder,subj.files[i]),
                                              eval            = eval, 
                                              ortho_detection = ortho_detection,
                                              comp_cores      = comp_cores, 
                                              ds.values       = ds.values,
                                              subject.id      = subject.id,
                                              ... )
             
                write.table(OrgQuery_vs_OrgSubj,
                            file.path(output.folder,paste0("map_q=",basename(query_file),"_s=",subj.files[i],"_orthodetec=",ortho_detection,"_eval=",eval,".csv")), 
                            sep       = sep,
                            col.names = TRUE,
                            row.names = FALSE,
                            quote     = FALSE)     
                
                if(progress.bar)
                        setTxtProgressBar(pb, i)
                
                }
        cat("\n")
        cat("\n")
        cat(paste0("All maps are stored in ",output.folder,"."))
        
}
        
        


