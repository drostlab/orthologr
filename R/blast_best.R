#' @title Perform a BLAST+ best hit search
#' @description This function performs a BLAST+ search (best hit) of a given set of protein sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}.
#' Default is \code{format} = \code{"fasta"}.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. 
#' \code{blast_algorithm} = \code{"blastp"}, \code{blast_algorithm} = \code{"blastn"}, \code{blast_algorithm} = \code{"tblastn"} .
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet. 
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' Please be aware that \code{max.target.seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param save.output a path to the location were the BLAST output shall be stored. E.g. \code{save.output} = \code{getwd()}
#' to store it in the current working directory, or \code{save.output} = \code{file.path(put,your,path,here)}.
#' @author Hajk-Georg Drost
#' @details Given a set of protein sequences (query sequences), a best hit blast search (BH BLAST) is being performed.
#' 
#' Internally to perform best hit searches, the BLAST+ parameter settings:
#' 
#' \code{"-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1"}
#' 
#' are used to speed up best hit computations.
#' 
#' @references
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schaeffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' @examples \dontrun{
#' 
#' # performing gene orthology inference using the best hit (BH) method
#' blast_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'            
#'            
#'            
#' # performing gene orthology inference using the best hit (BH) method starting with protein sequences
#' blast_best(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'            seq_type     = "protein")
#' 
#' 
#' 
#' # save the BLAST output file to the current working directory
#' blast_best(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'            seq_type     = "protein",
#'            save.output  = getwd())
#' 
#' 
#' 
#' 
#' 
#' # use multicore processing
#' blast_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores   = 2)
#'
#'
#'
#' # performing gene orthology inference using the best hit (BH) method and external
#' # blastp path
#' blast_best(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            path         = "path/to/blastp/")
#'            
#'            
#' }
#'
#' @return A data.table as returned by the \code{blast} function, storing the geneids
#' of orthologous genes (best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{blast}}, \code{\link{blast_rec}}, \code{\link{set_blast}}
#' @import data.table
#' @export
blast_best <- function(query_file, 
                       subject_file, 
                       seq_type        = "cds",
                       format          = "fasta", 
                       blast_algorithm = "blastp", 
                       delete_corrupt_cds = TRUE,
                       eval            = "1E-5",
                       max.target.seqs = 10000,
                       path            = NULL, 
                       comp_cores      = 1,
                       blast_params    = NULL, 
                       clean_folders   = FALSE,
                       save.output     = NULL){

        # default parameters for best hit filtering
         default_pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25"
        
        
        # performing a BLAST search from query against subject: blast(query,subject)
        
        hit_tbl <- blast( query_file      = query_file,
                             subject_file    = subject_file,
                             eval            = eval,
                             max.target.seqs = max.target.seqs,
                             delete_corrupt_cds = delete_corrupt_cds,
                             seq_type        = seq_type,
                             format          = format, 
                             path            = path, 
                             comp_cores      = comp_cores,
                             blast_params    = ifelse(!is.null(blast_params), paste0(blast_params," ",default_pars), default_pars),
                             clean_folders   = clean_folders,
                             save.output     = save.output)
        
        query_id <- '.' <- NULL
        hit_tbl <- dplyr::do(dplyr::group_by(hit_tbl, query_id), filter_best_hits(.))
        return(hit_tbl)
}

