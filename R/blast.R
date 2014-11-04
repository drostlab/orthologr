#' @title Interface function to BLAST+
#' @description This function performs a BLAST+ search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is \code{format} = \code{"fasta"}.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @details This function provides a fast communication between R and BLAST+. It is mainly used as internal functions
#' such as \code{\link{blast_best}} and \code{\link{blast_rec}} but can also be used to perform simple BLAST computations.
#' 
#' Note, that this function isn't as flexible as \code{\link{advanced_blast}}.
#' 
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @references 
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
#' \url{http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly}
#' 
#' \url{http://blast.ncbi.nlm.nih.gov/Blast.cgi}
#' @examples \dontrun{
#' # performing a BLAST search using blastp (default)
#' blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' 
#' # performing a BLAST search using blastp (default) using amino acid sequences as input file
#' blast(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type = "protein")
#'       
#' # in case you are working with a multicore machine, you can also run parallel
#' # BLAST computations using the comp_cores parameter: here with 2 cores
#' blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       comp_cores = 2)
#'       
#'
#'  # running blastp using additional parameters
#'  blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       blast_params = "-max_target_seqs 1")
#'
#'
#' # running blastp using additional parameters and an external blastp path
#'  blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       blast_params = "-max_target_seqs 1", path = "path/to/blastp/")              
#' }
#'
#' @return A data.table storing the BLAST hit table returned by BLAST.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
blast <- function(query_file, subject_file, seq_type = "cds",
                  format = "fasta", blast_algorithm = "blastp",
                  eval = "1E-5", path = NULL, comp_cores = 1,
                  blast_params = NULL, clean_folders = FALSE){
        
        if(!is.element(blast_algorithm,c("blastp")))
                stop("Please choose a valid BLAST mode.")
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        aa <- geneids <- NULL
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # initialize the BLAST search
        query.dt <- set_blast(file = query_file, seq_type = seq_type, format = format)[[1]]
        # make a BLASTable databse of the subject
        database <- set_blast(file = subject_file, seq_type = seq_type, format = format, makedb = TRUE)[[2]]
        # create an internal folder structure for the BLAST process 
        input = "blastinput.fasta"
        output = "blastresult.csv"
        
        
        if(!file.exists(paste0("_blast_db",f_sep))){
                
                dir.create("_blast_db")
        }
        
        currwd <- getwd()
        setwd(file.path(currwd, "_blast_db"))
        
        # determine the number of cores on a multicore machine
        cores <- parallel::detectCores()
        
        # in case one tries to use more cores than are available
        if(comp_cores > cores)
                stop("You chose more cores than are available on your machine.")
        
        write_AA <- as.list(query.dt[ ,aa])
        name <- query.dt[ ,geneids]      
        
        tryCatch(
         {
                seqinr::write.fasta(write_AA, names = name,
                                    nbchar = 80,open = "w",
                                    file.out = input)
         },error = function(){ stop(paste0("File ",input,
                                            " could not be written properly to the internal folder environment.",
                                            " Please check the path to ",input,".") ) }
         
        )
        
        # test whether the connection to BLAST+ works
        tryCatch(
        {      
                if(is.null(path)){
                        
                        if(blast_algorithm == "blastp"){
                                
                                if(is.null(blast_params)){
                                        
                                        # use the default parameters when running blastp
                                        system(
                                               paste0("blastp -db ",database," -query ",input,
                                                      " -evalue ",eval," -out ", output ," -outfmt 6", 
                                                      " -num_threads ", comp_cores)
                                               )
                                } else {
                                        
                                        # add additional parameters when running blastp
                                        system(
                                                paste0("blastp -db ",database," -query ",input,
                                                       " -evalue ",eval," -out ", output ," -outfmt 6", 
                                                       " -num_threads ", comp_cores," ",blast_params)
                                        )   
                                        
                                }
                        }
                        
                } else {
                        
                        if(blast_algorithm == "blastp"){
                                
                                if(is.null(blast_params)){
                                        
                                        # use the default parameters when running blastp
                                        system(
                                                paste0("export PATH=$PATH:",path,"; blastp -db ",
                                                       database," -query ",input," -evalue ",
                                                       eval," -out ", output ," -outfmt 6", 
                                                       " -num_threads ", comp_cores)
                                              )
                                } else {
                                        
                                        # add additional parameters when running blastp
                                        system(
                                                paste0("export PATH=$PATH:",path,"; blastp -db ",
                                                       database," -query ",input," -evalue ",
                                                       eval," -out ", output ," -outfmt 6", 
                                                       " -num_threads ", comp_cores," ",blast_params)
                                        )
                                        
                                        
                                }
                        }
                }
        
        },error = function(){ stop(paste0("Please check the correct path to ",blast_algorithm,
                                           "... the interface call did not work properly.") ) }
        )
        
        # additional blast parameters can be found here:
        # http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
        blast_table_names <- c("query_id","subject_id","perc_identity",
                               "alig_length","mismatches",
                               "gap_openings","q_start",
                               "q_end","s_start","s_end",
                               "evalue","bit_score")
        
        # define the colClasses for faster file streaming
        col_Classes <- c(rep("character",2), rep("numeric",10))
        
        tryCatch(
         {
                  hit_table <- data.table::fread(input = output, sep = "\t", 
                                                 header = FALSE, colClasses = col_Classes)
        
                  data.table::setnames(hit_table, old = paste0("V",1:length(blast_table_names)),
                                       new = blast_table_names)
        
                  data.table::setkeyv(hit_table, c("query_id","subject_id"))
                  
                  setwd(file.path(currwd))
                  
                  if(clean_folders)
                          clean_all_folders("_blast_db")
                  
                  return(hit_table)
         }, error = function(){ stop(paste0("File ",output, "could not be read correctly.",
                                             " Please check the correct path to ",output,
                                             " or whether BLAST did write the resulting hit table correctly.") ) }
        )

}











