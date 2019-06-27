#' @title Perform Divergence Stratigraphy
#' @description This function takes a query organism and performs
#' divergence stratigraphy (Quint et al.,2012 ; Drost et al. 2015) against a
#' closely related subject organism.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Available methods are: "BH" (BLAST best hit), "RBH" (BLAST reciprocal best hit).
#' @param dnds_est.method the dNdS estimation method that shall be used.
#' Options are:
#' \itemize{
#' \item \code{dnds_est.method = "Comeron"} (Default): Comeron's method (1995)
#' \item \code{dnds_est.method = "Li"}: Li's method (1993)
#' \item \code{dnds_est.method = "NG"}: Nei, M. and Gojobori, T. (1986)
#' \item \code{dnds_est.method = "LWL"}: Li, W.H., et al. (1985)
#' \item \code{dnds_est.method = "LPB"}: Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' \item \code{dnds_est.method = "MLWL"}: (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' \item \code{dnds_est.method = "YN"}: Yang, Z. and Nielsen, R. (2000)
#' \item \code{dnds_est.method ="MYN"} (Modified YN): Zhang, Z., et al. (2006)
#' }
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform
#'  parallel computations on a multicore machine.
#' @param dnds.threshold a numeric value specifying the dnds threshold for genes that shall be retained.
#' Hence all genes having a dNdS value <= \code{dnds.threshold} are retained. Default is \code{dnds.threshold} = 2.
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out to the console.
#' @param clean_folders a logical value specifying whether the internal folder structure shall be deleted (cleaned) after
#'  processing this function. Default is \code{clean_folders} = \code{FALSE}.
#' @param ds.values a logical value specifying whether divegrence stratum values (ds values) or dNdS values shall be returned
#' by \code{divergence_stratigraphy}. Default is \code{ds.values} = \code{TRUE}.
#' @param subject.id a logical value indicating whether \code{query_id} AND \code{subject_id} should be returned.
#' @details Introduced by Quint et al., 2012 and extended in Drost et al. 2015, divergence stratigraphy
#'  is the process of quantifying the selection pressure (in terms of amino acid sequence divergence) acting on
#'  orthologous genes between closely related species. The resulting sequence divergence map (short divergence map),
#'  stores the divergence stratum in the first column and the query_id of inferred orthologous genes in the second column.
#'  
#'  Following steps are performed to obtain a standard divergence map based on divergence_stratigraphy:
#'  
#'  Divergence Stratigraphy:
#'  
#'  \itemize{
#'  \item 1) Orthology Inference using BLAST reciprocal best hit ("RBH") based on blastp
#'  
#'  \item 2) Pairwise global amino acid alignments of orthologous genes using the Needleman-Wunsch algorithm
#'  
#'  \item 3) Codon alignments of orthologous genes using PAL2NAL
#'  
#'  \item 4) dNdS estimation using Comeron's method (1995)
#'  
#'  \item 5) Assigning dNdS values to divergence strata (deciles)
#' }
#' @note Although this function has been heavily optimized and parallelized, performing
#' Divergence Stratigraphy using two genomes will take some computation time.
#' 
#' In our experience performing Divergence Stratigraphy using two genomes (one query and one subject genome)
#' on an 8 core machine can take up to 1,5 - 2 hours.
#'   
#' @author Hajk-Georg Drost
#' @references
#'  
#'  Quint M et al. (2012). A transcriptomic hourglass in plant embryogenesis. Nature (490): 98-101.
#'  
#'  Drost HG et al. (2015). Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012
#'  
#' @examples \dontrun{
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'       
#'       
#'       
#'  # performing standard divergence stratigraphy using the blast_path argument to specify
#'  # the path to theblastp
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5", 
#'       ortho_detection = "RBH",
#'       blast_path      = "path/to/blastp",
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'  
#'  
#'  
#'  # Divergence Stratigraphy can also be performed in parallel 
#'  # (on a multicore machine) using the 'comp_cores' argument
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5", 
#'       ortho_detection = "RBH", 
#'       comp_cores      = 2, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE)
#'  
#'  
#'  
#'  
#'  # in case you want a divergence map with DS values and query_id and subject_id
#'  # you can specify subject.id = TRUE
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE,
#'       subject.id      = TRUE)
#'       
#'       
#'       
#'       
#'  # in case you want a divergence map with KaKs values and query_id and subject_id
#'  # you can specify ds.values = FALSE and subject.id = TRUE
#'  # performing standard divergence stratigraphy
#'  divergence_stratigraphy(
#'       query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       eval            = "1E-5",
#'       ortho_detection = "RBH", 
#'       comp_cores      = 1, 
#'       quiet           = TRUE, 
#'       clean_folders   = TRUE,
#'       subject.id      = TRUE)
#'  }
#'  
#' @return A data.table storing the divergence map of the query organism.
#' @seealso \code{\link{dNdS}}, \code{\link{divergence_map}},  \code{\link{substitutionrate}}, \code{\link{multi_aln}},
#'   \code{\link{codon_aln}}, \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{map.generator}}
#' @export
divergence_stratigraphy <- function(query_file, 
                                    subject_file, 
                                    eval            = "1E-5",
                                    ortho_detection = "RBH",
                                    dnds_est.method = "Comeron",
                                    delete_corrupt_cds = FALSE,
                                    blast_path      = NULL, 
                                    comp_cores      = 1,
                                    dnds.threshold  = 2,
                                    store_locally   = FALSE,
                                    quiet           = FALSE, 
                                    clean_folders   = FALSE, 
                                    ds.values       = TRUE,
                                    subject.id      = FALSE){
        
        if (!is.ortho_detection_method(ortho_detection))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        message("Running Divergence Stratigraphy ...")
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- dNdS <- NULL
        
        dNdS_tbl <- filter_dNdS( dNdS( query_file      = query_file,
                                       subject_file    = subject_file,
                                       ortho_detection = ortho_detection,
                                       delete_corrupt_cds = delete_corrupt_cds,
                                       aa_aln_type     = "pairwise", 
                                       aa_aln_tool     = "NW",
                                       codon_aln_tool  = "pal2nal", 
                                       dnds_est.method = dnds_est.method,
                                       eval            = eval,
                                       comp_cores      = comp_cores,
                                       store_locally   = store_locally,
                                       quiet           = quiet ), 
                                       dnds.threshold  = dnds.threshold)
        
        if (ds.values) {
                # divergence map: standard = col1: divergence stratum, col2: query_id
                dm_tbl <- divergence_map( dNdS_tbl = dNdS_tbl , subject.id = subject.id )
        }
        
        if (!ds.values) {
                if (!subject.id)
                        dm_tbl <- stats::na.omit(dNdS_tbl[ ,list(dNdS,query_id)]) 
                if (subject.id)
                        dm_tbl <- stats::na.omit(dNdS_tbl[ ,list(dNdS,query_id,subject_id)])
        }
        
        if (clean_folders)
                clean_all_folders(c(file.path(tempdir(),"_alignment"), file.path(tempdir(),"_blast_db"), file.path(tempdir(),"_calculation")))
        
        DivergenceMap.DF <- as.data.frame(dm_tbl)
        
        if (ds.values)
                colnames(DivergenceMap.DF)[1] <- "DS"
        if (!ds.values)
                colnames(DivergenceMap.DF)[1] <- "dNdS"
        
        message("Divergence Stratigraphy completed successfully.")
        return( DivergenceMap.DF )
}








