#' @title Compute dNdS values for two organisms
#' @description This function takes the CDS files of two organisms of interest (query_file and subject_file)
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms. 
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". See \code{\link{read.cds}},
#' \code{\link{read.genome}}, \code{\link{read.proteome}} for more details.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Further methods are: "BH" (BLAST best hit), "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), "OrthoMCL, "IP" (InParanoid).
#' @param cdd.path path to the cdd database folder (specify when using \code{ortho_detection} = \code{"DELTA"}).
#' @param blast_params a character string specifying additional parameters that shall be passed to BLAST. Default is \code{blast_params} = \code{NULL}. 
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param ortho_path a character string specifying the path to the orthology inference program such as \code{\link{ProteinOrtho}}, etc. (in case you don't use the default path).
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type} = "multiple" or \code{aa_aln_type} = "pairwise".
#' Default is \code{aa_aln_type} = "multiple".
#' @param aa_aln_tool a character string specifying the program that should be used e.g. "clustalw".
#' @param aa_aln_path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param aa_aln_params  a character string specifying additional parameters that shall be passed to the selected alignment tool. Default is \code{aa_aln_params} = \code{NULL} 
#' (no addintional parameters are passed to the selected alignment tool).
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' Right now only "pal2nal" can be selected as codon alignment tool.
#' @param kaks_calc_path a character string specifying the execution path to KaKs_Calculator. Default is \code{kaks_calc_path} = \code{NULL}
#' (meaning that KaKs_Calculator is stored and executable in your default \code{PATH}).
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine. 
#' @param quiet a logical value specifying whether the output of the corresponding alignment tool shall be printed out to the console.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @author Hajk-Georg Drost
#' @details 
#' 
#' 
#' The dN/dS ratio quantifies the mode and strength of selection acting on a pair of orthologous genes.
#' This selection pressure can be quantified by comparing synonymous substitution rates (dS) that are assumed to be neutral
#' with nonsynonymous substitution rates (dN), which are exposed to selection as they 
#' change the amino acid composition of a protein (Mugal et al., 2013 http://mbe.oxfordjournals.org/content/31/1/212).
#' 
#' The \pkg{orthologr} package provides the \code{\link{dNdS}} function to perform dNdS estimation on pairs of orthologous genes.
#' This function takes the CDS files of two organisms of interest (\code{query_file} and \code{subject_file}) 
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms.
#' 
#' The following pipieline resembles the dNdS estimation process:
#' 
#' \itemize{
#'                 
#' \item 1) Orthology Inference: e.g. BLAST reciprocal best hit (RBH)
#' 
#' \item 2) Pairwise sequence alignment: e.g. clustalw for pairwise amino acid sequence alignments
#' 
#' \item 3) Codon Alignment: e.g. pal2nal program
#' 
#' \item 4) dNdS estimation: e.g. Yang, Z. and Nielsen, R. (2000) \url{http://mbe.oxfordjournals.org/content/17/1/32.short}
#' 
#' }
#' Note: it is assumed that when using \code{dNdS()} all corresponding multiple sequence alignment programs you
#' want to use are already installed on your machine and are executable via either
#' the default execution \code{PATH} or you specifically define the location of the executable file
#' via the \code{aa_aln_path} or \code{blast_path} argument that can be passed to \code{dNdS()}.
#' 
#' The \code{dNdS()} function can be used choosing the folllowing options:
#' 
#' \itemize{
#' \item \code{ortho_detection} : 
#'    \itemize{  
#'    \item \code{"RBH"} (BLAST best reciprocal hit)
#'    \item \code{"BH"} (BLAST best hit)
#'    \item \code{"DELTA"} (DELTA-BLAST best reciprocal hit) 
#'    \item \code{"PO"} (ProteinOrtho) 
#'      }
#'      
#' \item \code{aa_aln_type} : 
#'  \itemize{
#'   \item \code{"multiple"}
#'   \item \code{"pairwise"}
#'  }
#' 
#' \item \code{aa_aln_tool} :
#' \itemize{
#'  \item \code{"clustalw"}
#'  \item \code{"t_coffee"}
#'  \item \code{"muscle"}
#'  \item \code{"clustalo"}
#'  \item \code{"mafft"}
#'  \item \code{"NW"} (in case \code{aa_aln_type = "pairwise"})
#' }
#' 
#' \item \code{codon_aln_tool} :
#' \itemize{
#'  \item \code{"pal2nal"}
#'  }
#' \item \code{dnds_est.method} : 
#' \itemize{
#' \item "Li" : Li's method (1993)
#' \item "Comeron" : Comeron's method (1995)
#' \item "NG": Nei, M. and Gojobori, T. (1986)
#' \item "LWL": Li, W.H., et al. (1985)
#' \item "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' \item "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' \item "YN": Yang, Z. and Nielsen, R. (2000)
#' \item "MYN" (Modified YN): Zhang, Z., et al. (2006)
#' 
#' }
#' 
#' }
#' 
#' @references
#' 
#' seqinr: \url{http://seqinr.r-forge.r-project.org/}
#' 
#' Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator: Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263. 
#' 
#' KaKs_Calculator: \url{https://code.google.com/p/kaks-calculator/} [GNU GPL-3 license]
#' 
#' Paradis, E. (2012) Analysis of Phylogenetics and Evolution with R (Second Edition). New York: Springer.
#'
#' Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289-290.
#' 
#' More information on \pkg{ape} can be found at \url{http://ape-package.ird.fr/}.
#' 
#' Pages H, Aboyoun P, Gentleman R and DebRoy S. Biostrings: String objects representing biological sequences, and
#' matching algorithms. R package version 2.32.1.
#' 
#' @return A data.table storing the dNdS values of the correspnding genes.
#' @examples \dontrun{
#' 
#' # get a dNdS table using:
#' # 1) reciprocal best hit for orthology inference (RBH)
#' # 2) Needleman-Wunsch for pairwise amino acid alignments
#' # 3) pal2nal for codon alignments
#' # 4) Comeron for dNdS estimation
#' # 5) single core processing 'comp_cores = 1'
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH", 
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW", 
#'      codon_aln_tool  = "pal2nal", 
#'      dnds_est.method = "Comeron", 
#'      comp_cores      = 1 )
#' 
#' 
#' # running dNdS using the 'aa_aln_path' argument to specify the path to
#' # the corresponding alignment tool
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH", 
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW", 
#'      aa_aln_path     = "/path/to/clustalw/",
#'      codon_aln_tool  = "pal2nal", 
#'      dnds_est.method = "Comeron", 
#'      comp_cores      = 1 )
#' 
#' # The same result can be obtained using multicore processing using: comp_cores = 2 or 3 or more ...
#' dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'      ortho_detection = "RBH", 
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW", 
#'      aa_aln_path     = "/path/to/clustalw/",
#'      codon_aln_tool  = "pal2nal", 
#'      dnds_est.method = "Comeron", 
#'      comp_cores      = 1 )
#' 
#' }
#' @seealso \code{\link{divergence_stratigraphy}}, \code{\link{orthologs}}, \code{\link{substitutionrate}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, \code{\link{blast_best}},
#' \code{\link{blast_rec}}, \code{\link{read.cds}} 
#' @export

dNdS <- function(query_file, 
                 subject_file, 
                 seq_type        = "protein",
                 format          = "fasta", 
                 ortho_detection = "RBH",
                 cdd.path        = NULL,
                 blast_params    = NULL, 
                 blast_path      = NULL, 
                 eval            = "1E-5", 
                 ortho_path      = NULL, 
                 aa_aln_type     = "multiple", 
                 aa_aln_tool     = "clustalw", 
                 aa_aln_path     = NULL, 
                 aa_aln_params   = NULL, 
                 codon_aln_tool  = "pal2nal", 
                 kaks_calc_path  = NULL, 
                 dnds_est.method = "YN", 
                 comp_cores      = 1, 
                 quiet           = FALSE, 
                 clean_folders   = FALSE){
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if (!is.ortho_detection_method(ortho_detection))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        if (!is.element(aa_aln_type, c("multiple", "pairwise")))
                stop("Please choose a supported alignement type: 'multiple' or 'pairwise'")
        
        if (aa_aln_type == "multiple") {
                if (!is.multiple_aln_tool(aa_aln_tool))
                        stop(
                                "Please choose a multiple alignment tool that is supported by this function or try to choose aa_aln_type = 'pairwise'."
                        )
        }
        
        if (aa_aln_type == "pairwise") {
                if (!is.pairwise_aln_tool(aa_aln_tool))
                        stop(
                                "Please choose a pairwise alignment tool that is supported by this function or try to choose aa_aln_type = 'multiple'."
                        )
        }
        
        if (!is.codon_aln_tool(codon_aln_tool))
                stop("Please choose a codon alignment tool that is supported by this function.")
        
        if (!is.dnds_est_method(dnds_est.method))
                stop("Please choose a dNdS estimation method that is supported by this function.")
        
        aa <- geneids <- NULL
        
        # blast each translated aminoacid sequence against the related database to get a 
        # hit table with pairs of geneids  
        
        # use BLAST best hit as orthology inference method
        if (ortho_detection == "BH") {
                # seq_type = "cds" -> dNdS() needs CDS files as input!
                hit.table <- data.table::copy(
                        blast_best(
                                query_file   = query_file,
                                subject_file = subject_file,
                                blast_params = blast_params,
                                path         = blast_path,
                                comp_cores   = comp_cores,
                                seq_type     = "cds",
                                eval         = eval,
                                format       = format
                        )
                )
                
                q_cds <- read.cds(file   = query_file,
                                  format = format)
                
                s_cds <- read.cds(file   = subject_file,
                                  format = format)
                
                filename_qry <-
                        unlist(strsplit(
                                query_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_qry <- filename_qry[length(filename_qry)]
                
                input = paste0("query_", filename_qry, ".fasta")
                
                q_aa <-
                        read.proteome(file = file.path(tempdir(), "_blast_db", input),
                                      format = "fasta")
                
                s_aa_tmp <- cds2aa(s_cds)
                
                filename_subj <-
                        unlist(strsplit(
                                subject_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_subj <-
                        filename_subj[length(filename_subj)]
                
                seqinr::write.fasta(
                        sequences = as.list(s_aa_tmp[, aa]),
                        names     = s_aa_tmp[, geneids],
                        nbchar    = 80,
                        open      = "w",
                        file.out  = file.path(tempdir(), "_blast_db", filename_subj)
                )
                
                s_aa <-
                        read.proteome(file = file.path(tempdir(), "_blast_db", filename_subj),
                                      format = "fasta")
                
        }
        
        # use BLAST best reciprocal hit as orthology inference method
        if (ortho_detection == "RBH") {
                # seq_type = "cds" -> dNdS() needs CDS files as input!
                hit.table <- data.table::copy(
                        blast_rec(
                                query_file   = query_file,
                                subject_file = subject_file,
                                blast_params = blast_params,
                                path         = blast_path,
                                comp_cores   = comp_cores,
                                seq_type     = "cds",
                                eval         = eval,
                                format       = format
                        )
                )
                
                
                q_cds <- read.cds(file   = query_file,
                                  format = format)
                
                s_cds <- read.cds(file   = subject_file,
                                  format = format)
                
                
                filename_qry <-
                        unlist(strsplit(
                                query_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_qry <- filename_qry[length(filename_qry)]
                
                input_qry = paste0("query_", filename_qry, ".fasta")
                
                q_aa <-
                        read.proteome(file = file.path(tempdir(), "_blast_db", input_qry),
                                      format = "fasta")
                
                filename_subj <-
                        unlist(strsplit(
                                subject_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_subj <-
                        filename_subj[length(filename_subj)]
                
                input_subj = paste0("query_", filename_subj, ".fasta")
                
                s_aa <-
                        read.proteome(file = file.path(tempdir(), "_blast_db", input_subj),
                                      format = "fasta")
        
        }
        
        if (!is.element(ortho_detection, c("BH", "RBH"))) {
                # seq_type = "cds" -> dNdS() needs CDS files as input!
                hit.table <- data.table::copy(
                        orthologs(
                                query_file      = query_file,
                                subject_files   = subject_file,
                                ortho_detection = ortho_detection,
                                eval            = eval,
                                path            = ortho_path,
                                cdd.path        = cdd.path,
                                comp_cores      = comp_cores,
                                seq_type        = "cds",
                                format          = format,
                                detailed_output = FALSE,
                                quiet           = quiet,
                                clean_folders   = FALSE
                        )
                )
                
                
                q_cds <- read.cds(file   = query_file,
                                  format = format)
                
                s_cds <- read.cds(file   = subject_file,
                                  format = format)
                
                
                q_aa_tmp <- cds2aa(q_cds)
                s_aa_tmp <- cds2aa(s_cds)
                
                filename_qry <-
                        unlist(strsplit(
                                query_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_qry <- filename_qry[length(filename_qry)]
                
                filename_subj <-
                        unlist(strsplit(
                                subject_file,
                                f_sep,
                                fixed = FALSE,
                                perl = TRUE,
                                useBytes = FALSE
                        ))
                filename_subj <-
                        filename_subj[length(filename_subj)]
                
                seqinr::write.fasta(
                        sequences = as.list(q_aa_tmp[, aa]),
                        names     = q_aa_tmp[, geneids],
                        nbchar    = 80,
                        open      = "w",
                        file.out  = file.path(tempdir(), filename_qry)
                )
                
                seqinr::write.fasta(
                        sequences = as.list(s_aa_tmp[, aa]),
                        names     = s_aa_tmp[, geneids],
                        nbchar    = 80,
                        open      = "w",
                        file.out  = file.path(tempdir(), filename_subj)
                )
                
                q_aa <-
                        read.proteome(file = file.path(tempdir(), filename_qry),
                                      format = "fasta")
                s_aa <-
                        read.proteome(file = file.path(tempdir(), filename_subj),
                                      format = "fasta")
                
        }
        
        
        data.table::setnames(q_cds, old = c("geneids", "seqs"), new = c("query_id","query_cds"))
        data.table::setnames(s_cds, old = c("geneids", "seqs"), new = c("subject_id","subject_cds"))
        data.table::setnames(q_aa, old = c("geneids", "seqs"), new = c("query_id","query_aa"))
        data.table::setnames(s_aa, old = c("geneids", "seqs"), new = c("subject_id","subject_aa"))
        
        # joining all tables to a final table containing: query_id, subject_id, query_aa_seq, subject_aa_seq, query_cds_seq, and subject_cds_seq
        query_tbl <- dplyr::inner_join(dtplyr::tbl_dt(q_cds), dtplyr::tbl_dt(q_aa), by = "query_id")
        subject_tbl <- dplyr::inner_join(dtplyr::tbl_dt(s_cds), dtplyr::tbl_dt(s_aa), by = "subject_id")
        joint_query_tbl <- dplyr::inner_join(dtplyr::tbl_dt(hit.table), dtplyr::tbl_dt(query_tbl), by = "query_id")
        joint_subject_tbl <- dplyr::inner_join(dtplyr::tbl_dt(hit.table), dtplyr::tbl_dt(subject_tbl), by = "subject_id")
        final_tbl <- dplyr::inner_join(dtplyr::tbl_dt(joint_query_tbl), dtplyr::tbl_dt(joint_subject_tbl), by = c("query_id","subject_id"))
        
        if (nrow(final_tbl) == 0)
                stop("No orthologs could be found! Please check your input files!")
        
        
       # return the dNdS table for all query_ids and subject_ids
           return( compute_dnds( complete_tbl    = final_tbl,
                                 aa_aln_type     = aa_aln_type,
                                 aa_aln_tool     = aa_aln_tool,
                                 aa_aln_path     = aa_aln_path,
                                 codon_aln_tool  = codon_aln_tool, 
                                 kaks_calc_path  = kaks_calc_path, 
                                 dnds_est.method = dnds_est.method, 
                                 quiet           = quiet,
                                 comp_cores      = comp_cores, 
                                 clean_folders   = clean_folders )
                  )   
}






