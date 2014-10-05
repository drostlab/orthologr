#' @title Interface function to BLAST+
#' @description This function performs a BLAST+ search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @references 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#' 
#' http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi
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

#'  # running blastp using additional parameters
#'  blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       blast_params = "-max_target_seqs 1")
#'              
#' }
#'
#' @return A data.table storing the BLAST hit table returned by BLAST.
#' @import seqinr
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
blast <- function(query_file, subject_file, seq_type = "cds",
                  format = "fasta", blast_algorithm = "blastp",
                  eval = "1E-5", path = NULL, comp_cores = 1,
                  blast_params = NULL){
        
        if(!is.element(blast_algorithm,c("blastp")))
                stop("Please choose a valid BLAST mode.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # initialize the BLAST search
        query.dt <- set_blast(file = query_file, seq_type = seq_type, format = format)[[1]]
        # make a BLASTable databse of the subject
        database <- set_blast(file = subject_file, seq_type = seq_type, format = format, makedb = TRUE)[[2]]
        # create an internal folder structure for the BLAST process 
        input = paste0("_blast_db",f_sep,"blastinput.fasta") 
        output = paste0("_blast_db",f_sep,"blastresult.csv")
        
        
        if(!file.exists(paste0("_blast_db",f_sep))){
                
                dir.create("_blast_db")
        }
        
        
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
        
        },error = function(){ stop(paste0("Please check the correct path to ",tool,
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
                  
                  
                  return(hit_table)
         }, error = function(){ stop(paste0("File ",output, "could not be read correctly.",
                                             " Please check the correct path to ",output,
                                             " or whether BLAST did write the resulting hit table correctly.") ) }
        )
}


#' @title Function to perform a BLAST+ best hit search
#' @description This function performs a blast+ search (best hit) of a given set of protein sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details Given a set of protein sequences A, a best hit blast search is being performed from A to database.
#' @references
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' @examples \dontrun{
#' 
#' # performing gene orthology inference using the best hit (BH) method
#' blast_best(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'            
#' # performing gene orthology inference using the best hit (BH) method starting with protein sequences
#' blast_best(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'            seq_type = "protein")
#' 
#' # use multicore processing
#' blast_best(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores = 2)
#'            
#' }
#'
#' @return A data.table as returned by the \code{blast} function, storing the geneids
#' of orthologous genes (best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{blast}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
blast_best <- function(query_file, subject_file, seq_type = "cds",
                       format = "fasta", blast_algorithm = "blastp", 
                       eval = "1E-5", path = NULL, comp_cores = 1,
                       blast_params = NULL){
        
        # default parameters for best hit filtering
        default_pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1"
        
        
        # performing a BLAST search from query against subject: blast(query,subject)
        # using the BLAST parameter: '-max_target_seqs 1' allows to retain
        # only the best hit result and therefore, speeds up the BLAST search process
        hit_tbl.dt <- blast(query_file = query_file, 
                            subject_file = subject_file,
                            seq_type = seq_type,
                            format = format, path = path, 
                            comp_cores = comp_cores,blast_params = ifelse(!is.null(blast_params),
                                                                          paste0(blast_params,
                                                                          " ",default_pars),
                                                                          default_pars))
       
        tryCatch(
         {
                 besthit_tbl <- hit_tbl.dt[ , sapply(.SD[ , evalue],min)[1],
                                by = key(hit_tbl.dt)]
        
                 #data.table::setnames(besthit_tbl, old = c("V1","V2"), new = c("subject_id","evalue"))
                 data.table::setnames(besthit_tbl, old = c("V1"), new = c("evalue"))
       
                 data.table::setkeyv(besthit_tbl, c("query_id","subject_id"))
                 
                 
                 # return a data.table storing only the best hits from the resulting 
                 # BLAST search
                 return( besthit_tbl )
         }, error = function() {stop(paste0("The BLAST output couldn't be read properly, maybe a problem occured when 
                                       selecting best hits from the resulting BLAST hit table."))} 
        )
              
}


#' @title Function to perform a BLAST+ reciprocal best hit (RBH) search
#' @description This function performs a blast+ search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details Given a set of protein sequences A and a different set of protein sequences B,
#' first a best hit blast search is being performed from A to B: blast(A,B) and afterwards
#' a best hit blast search is being performed from B to A: blast(B,A). Only protein sequences
#' that were found to be best hits in both directions are retained and returned.
#'
#' @references
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' @examples \dontrun{
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' blast_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#'           
#' # performing gene orthology inference using the reciprocal best hit (RBH) method
#' # starting with protein sequences
#' blast_rec(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'           seq_type = "protein")
#' 
#' # use multicore processing
#' blast_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores = 2)
#' }
#'
#' @return A data.table as returned by the blast() function, storing the geneids
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{blast}}, \code{\link{blast_best}}, \code{\link{advanced_blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
blast_rec <- function(query_file, subject_file, seq_type = "cds",
                      format = "fasta", blast_algorithm = "blastp", 
                      eval = "1E-5", path = NULL, comp_cores = 1, 
                      blast_params = NULL){
        
        orthoA <- blast_best(query_file,subject_file, 
                             format = format, seq_type = seq_type,
                             blast_algorithm = blast_algorithm,
                             path = path, comp_cores = comp_cores, blast_params = blast_params)
        
        orthoB <- blast_best(subject_file,query_file,
                             seq_type = seq_type,
                             format = format,blast_algorithm = blast_algorithm,
                             path = path, comp_cores = comp_cores, blast_params = blast_params)
        
        data.table::setnames(orthoB, old = c("query_id","subject_id"), new = c("subject_id","query_id"))
        
        tryCatch(
        {       
                
                return ( dplyr::semi_join(dplyr::tbl_dt(orthoA), dplyr::tbl_dt(orthoB),
                                          by = c("query_id","subject_id")) )
                
        }, error = function(){ stop(paste0("The BLAST tables resulting from ",query_file, " and ",
        subject_file," could not be joined properly to select only the reciprocal best hits."))}
        )
}


#' @title Function preparing the parameters and databases for subsequent BLAST searches
#' @description This function reads a file storing a specific sequence type, such as "cds", "protein", or
#' "dna" in a standard sequence file format such as "fasta", etc. and depending of the \code{makedb}
#'  parameter either creates a blast-able database, or returns the corresponding protein sequences
#'  as data.table object for further BLAST searches.  
#' @param file a character string specifying the path to the file storing the sequences of interest.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param makedb TRUE or FALSE whether a database should be created or not (BLAST parameter 'makeblastdb').
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param makedb_type a character string specifying the sequence type stored in the BLAST database
#' that is generated using 'makeblastdb'. Options are: "protein" and "nucleotide". Default is \code{makedb_type} = "protein".
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @return A list storing two elements. The first element [[1]] corresponds to the data.table storing the gene ids in the first column and
#' the corresponding dna (cds) sequence in the second column and the aminoacid sequence third column.
#' The second list element [[2]] stores the name of the protein database that was created by 'makeblastdb'.
#' @references
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#'
#' @examples \dontrun{
#'  # running the set function to see an example output
#'  head(set_blast(file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'))[[1]] , 2)
#' }
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{blast}}, \code{\link{advanced_makedb}}
#' @export
set_blast <- function(file, seq_type = "cds",format = "fasta", makedb = FALSE,
                      path = NULL, makedb_type = "protein", ...){
        
        # HERE WE NEED TO INSERT SOME QUALITY CONTROL
        # CONCERNING THE CASE THAT A POTENTIAL USER
        # COULD INSERT SOME NON-CDS SEQUENCES
        
        if(!is.element(seq_type,c("cds","protein","dna")))
                stop("Please choose either: 'cds', 'protein', or 'dna' as seq_type.")
        
        if(!is.element(makedb_type,c("protein","nucleotide")))
                stop("Please choose either: 'protein' or 'nucleotide' as BLAST database type (makedb_type).")
        
        if(makedb_type == "protein")
                db_type <- "prot"
        
        if(makedb_type == "nucleotide")
                db_type <- "nucl"
        
        if(seq_type == "cds"){
                # read cds file
                # copy the data.table due to this discussion:
                # http://stackoverflow.com/questions/8030452/pass-by-reference-the-operator-in-the-data-table-package
                dt <- data.table::copy(read.cds(file = file, format = "fasta", ...))
        
                if(!is.data.table(dt))
                        stop("Your CDS file was not corretly transformed into a data.table object.")
                
                # When using data.tables within packaes, always make sure
                # 'data.table' is included in the DESCRIPTION file as 'Imports' AND
                # in the NAMESPACE file with 'imports(data.table)' -> @import when using roxygen2
                # http://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package
                # https://github.com/hadley/dplyr/issues/548
                
                # omit empty sequences
                dt <- dt[,.SD[sapply(seqs,function(x){return(! (is.na(x) || x=="") )})]]
                
                # omit sequences taht are not multiples of 3
                dt <- dt[,.SD[sapply(seqs,function(x){return(nchar(x)%%3==0)})]]

                # omit sequences consisting of others than ACGT
                dt <- dt[,.SD[sapply(seqs,is.dnaSequence)]]
                
                # translate cds to protein sequences
                tryCatch(
                         {
                                 dt[ , aa := transl(seqs), by = geneids]
        
                          }, error = function() {stop(paste0("The input coding sequences could not be translated properly to amino acid sequences.",
                                   ,"\n"," Please check whether ",file, " stores valid coding sequences."))}
                          )

        }
        
        
        # read file storing protein sequence -> no translation needed
        if(seq_type == "protein"){
                
                dt <- data.table::copy(read.proteome(file = file, format = "fasta", ...))
                data.table::setnames(dt, old = "seqs", new= "aa")
        }
        
        
        if(seq_type == "dna"){
                
                dt <- data.table::copy(read.genome(file = file, format = "fasta", ...))
                
                # here will come -> CDS prediction
                # then CDS -> protein translation
                
        }
        
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        
        # WE COULD ALSO INSERT A IF-STATEMENT CHECKING THAT
        # IN CASE THE USER INSERTS AN AA-FASTA FILE,
        # THE TRANSLATION PROCESS IS BEING OMITTED
        
        
        # makedb
        dbname <- vector(mode = "character", length = 1)
        filename <- unlist(strsplit(file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
        filename <- filename[length(filename)]
        
        
        if(makedb){
                if(!file.exists(paste0("_blast_db",f_sep))){ 
                        
                        dir.create("_blast_db") 
                
                }
                
                dbname <- paste0("_blast_db",f_sep,"out_",filename,"_translate.fasta")
                seqinr::write.fasta(as.list(dt[ , aa]),
                                    names = dt[ , geneids],
                                    nbchar = 80, open = "w",
                                    file.out = dbname 
                                    )
                
                tryCatch(
                {
                        if(is.null(path)){
                                system(
                                        paste0("makeblastdb -in ", dbname,
                                               " -input_type fasta -dbtype ",db_type," -hash_index")
                                       )
                        } else {
                                system(
                                        paste0("export PATH=",path,"; makeblastdb -in ",
                                               dbname," -input_type fasta -dbtype ",db_type," -hash_index")
                                      )
                        }
                }, error = function(){ stop(paste0("makeblastdb did not work properly. The default parameters are: \n",
                "-input_type fasta -dbtype prot . \n","Please check that you really want to work with a protein database.\n",
                "Additionally check: ",dbname," ."))}
                )
        }
        
        return(list(na.omit(dt), dbname))
}



#' @title Advanced interface function to BLAST+
#' @description This function performs a BLAST+ search of a given set of sequences and a given set of parameters against a BLAST database.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn","deltablast",... .
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param makedb_type a character string specifying the sequence type stored in the BLAST database
#' that is generated using 'makeblastdb'. Options are: "protein" and "nucleotide". Default is \code{makedb_type} = "protein".
#' @param taxonomy a logical value specifying whether the subject taxonomy shall be returned based on NCBI Taxonomy queries. Default is \code{FALSE}.
#' @param db_path a character string specidying the path to the local BLAST database.
#' @param sql_database a logical value specifying whether an SQL database shall be created to store the BLAST output.
#' This is only useful when the amount of data does not fit in-memory anymore.
#' @param session_name a character string specifying the name of the BLAST output. Default is \code{session_name} = \code{NULL}.
#' @param write.only a logical value specifying whether the resulting blast output should be returned as
#' data.table or database connection. Default is \code{write.only} = \code{FALSE}.
#' @details
#'  
#' Following BLAST programs and algorithms can be assigned to \code{blast_algorithm}:
#' 
#' "blastp", "blastn", "megablast","psi-blast", "phi-blast", "delta-blast", "blastx", "tblastn", "tblastx"
#' 
#' 
#' The intention of this function is to provide the user with an easy to use interface function
#' to BLAST which uses the same notation as the original BLAST program.
#' 
#' When using \code{taxonomy} =  \code{TRUE}, then the taxdb.btd and taxdb.bti must be stored
#' in the same folder as BLAST-able database. In case you need to specify the path to your BLAST database,
#' use the \code{db_path} argument.
#' 
#' The \code{advanced_blast} function can also store the BLAST output CSV file in an SQLite database.
#' This works only with dplyr version >= 0.3 . To store the BLAST output in an SQLite database and
#' to receive an SQLite connection as specified by \code{tbl} in \code{dplyr} please use the \code{sql_database} = \code{TRUE} argument. 
#' 
#' 
#' @author Hajk-Georg Drost
#' @note This function is also able to return the subject taxonomy.
#' 
#' For this purpose the NCBI taxdb must be installed on your system:
#' 
#' ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
#' 
#' In case you want to use the taxonomy option, set \code{taxonomy} = \code{TRUE}.
#' @references 
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi
#' @examples \dontrun{
#' 
#' 
#' # performing a BLAST search using blastp (default)
#' advanced_blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'       blast_algorithm = "blastp", blast_params = "-evalue 1E-5 -num_threads 2")
#'       
#' # performing a BLAST search using blastp (default) and starting with protein sequences
#' advanced_blast(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'       seq_type = "protein", blast_algorithm = "blastp", blast_params = "-evalue 1E-5 -num_threads 2")       
#'  
#'              
#' # when performing an advanced BLAST search, you can easily select the best hit using
#' library(dplyr)
#' 
#' adv_blast_test <- advanced_blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                                  subject_file = system.file('seqs/ortho_lyra_cds_1000.fasta', package = 'orthologr'),
#'                                  seq_type = "cds",blast_algorithm = "blastp", blast_params = "-evalue 1E-5 -num_threads 1")
#'                                  
#' best_hit <- adv_blast_test %>% group_by(query_id) %>% summarise(min(evalue))
#'
#'
#' # using a SQLite database to store the BLAST output and
#' # select some example data
#' 
#' library(dplyr)
#' 
#' sql_example <- advanced_blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                               subject_file = system.file('seqs/ortho_lyra_cds_1000.fasta', package = 'orthologr'),
#'                               seq_type = "cds",blast_algorithm = "blastp", blast_params = "-evalue 1E-5 -num_threads 1",
#'                               sql_database = TRUE)
#'
#' head(sql_example)
#' 
#' glimpse(sql_example)
#' 
#' # select all rows that have an evalue of zero
#' filter(sql_example,evalue == 0)
#' 
#' # select the best hit using the evalue criterion
#' sql_example %>% group_by(query_id) %>% summarise(best_hit_eval = min(evalue))
#'                                                  
#' }
#'
#' @return the \code{advanced_blast} function creates a folder named '_blast_db' and stores
#' the hit table returned by BLAST in this folder.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
advanced_blast <- function(query_file, subject_file, 
                           seq_type = "cds",format = "fasta", 
                           blast_algorithm = "blastp", path = NULL,
                           blast_params = NULL, makedb_type = "protein",
                           taxonomy = FALSE, db_path = NULL,
                           sql_database = FALSE, session_name = NULL,
                           write.only = FALSE){
        
        # http://blast.ncbi.nlm.nih.gov/Blast.cgi
        if(!is.element(blast_algorithm,c("blastp","blastn", "megablast",
                                         "psiblast", "phiblast", "deltablast",
                                         "blastx","tblastn","tblastx")))
                stop("Please choose a valid BLAST mode.")
        
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # initialize the BLAST search
        
        query.dt <- set_blast(file = query_file, format = format, 
                              seq_type = seq_type)[[1]]
        # make a BLASTable databse of the subject
        database <- set_blast(file = subject_file, format = format, seq_type = seq_type, 
                              makedb = TRUE, makedb_type = makedb_type)[[2]]
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(is.null(session_name))
                session_name <- get_filename(query_file)
        
        # create an internal folder structure for the BLAST process 
        input = paste0("_blast_db",f_sep,session_name,"_blastinput.fasta") 
        output = paste0("_blast_db",f_sep,session_name,"_blastresult.csv")
        
        
        if(!file.exists(paste0("_blast_db",f_sep))){
                
                dir.create("_blast_db")
        }
        
        # write query fasta file 
        write_AA <- as.list(query.dt[ ,aa])
        name <- query.dt[ , geneids]      
        
        tryCatch(
         {
                 seqinr::write.fasta(write_AA, names = name,
                                     nbchar = 80,open = "w",
                                     file.out = input)
                 
         }, error = function(){ stop(paste0("File ",input," could not be written properly to the internal folder environment.\n",
         "Please check: ",query_file, " and ",subject_file,"."))}
        )
        
        if(is.null(path)){
                
                # here: http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
                # in column -outfmt additional output columns can be selected: ' 6 qseqid sseqid staxids sskingdoms' .. or ' 6 std'
                # which is the default value and the same as: ' 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
                tryCatch(
                               {
                                       if(taxonomy == TRUE){
                                
                                               if(!is.null(db_path)){
                                                       system( paste0("export BLASTDB=",db_path,"; ",blast_algorithm," -db ",database," -query ",input,
                                                               " -out ", output ," ",blast_params,
                                                               " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                               length mismatch gapopen qstart qend sstart send evalue 
                                                               bitscore score qcovs'")
                                                              )
                                                }
                                               
                                               if(is.null(db_path)){
                                                       system( paste0(blast_algorithm," -db ",database," -query ",input,
                                                                      " -out ", output ," ",blast_params,
                                                                      " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                               length mismatch gapopen qstart qend sstart send evalue 
                                                               bitscore score qcovs'")
                                                       )
                                               }
                                               
                                        }
                              }, error = function(){ stop(paste0("taxdb could not be included to the BLAST search. \n",
                                                     "Please check the validity of the path: ",db_path," .\n",
                                                     "Additionally, check the validity of: ",database,", ",input,
                                                     ", ",output,", and blast_params: ",blast_params," ."))}
                      )
                      
                tryCatch(
                        {      
                                if(taxonomy == FALSE){
                                
                                         system( paste0(blast_algorithm," -db ",database," -query ",input,
                                                 " -out ", output ," ",blast_params, 
                                                 " -outfmt '6 qseqid sseqid pident nident length mismatch 
                                                 gapopen qstart qend sstart send evalue bitscore score qcovs'")
                                               )
                                
                                }
                        }, error = function(){stop(paste0("The advanced BLAST search did not work properly.\n",
                        "Please check the validity of: ",database,", ",input,", ",output,", and blast_params: ",blast_params," ."))}
                )
                
        } else {
                
                tryCatch(
                         {
                                 if(taxonomy == TRUE){
                        
                                         if(!is.null(db_path)){     
                                
                                                 system(
                                                        paste0("export PATH=",path,"; ","export BLASTDB=",db_path,"; ",blast_algorithm," -db ",
                                                        database," -query ",input," -out ", output ," ",blast_params,
                                                        " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                        length mismatch gapopen qstart qend sstart send evalue bitscore 
                                                        score qcovs'")
                                                       )
                                         } 
                                         
                                         
                                         if(is.null(db_path)){
                                                 
                                                 system(
                                                         paste0("export PATH=",path,"; ",blast_algorithm," -db ",
                                                                database," -query ",input," -out ", output ," ",blast_params,
                                                                " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                        length mismatch gapopen qstart qend sstart send evalue bitscore 
                                                        score qcovs'")
                                                 )
                                                 
                                                 
                                         }
                       
                                  }
                         }, error = function(){stop(paste0("taxdb could not be included to the BLAST search. \n",
                                                            "Please check the validity of the path: ",db_path," .\n",
                                                            "Additionally, check the validity of: ",database,", ",input,",
                                                            ",output,", and blast_params: ",blast_params," ."))}
                )
                
                tryCatch(
                         {
                                 if(taxonomy == FALSE){
                        
                                         system(
                                                paste0("export PATH=$PATH:",path,"; ",blast_algorithm," -db ",
                                                database," -query ",input," -out ", output ," ",blast_params,
                                                " -outfmt '6 qseqid sseqid pident nident length 
                                                mismatch gapopen qstart qend sstart send evalue bitscore 
                                                score qcovs'")
                                                )
                        
                                 }
                         }, error = function(){stop(paste0("The advanced BLAST search did not work properly.\n",
                                                            "Please check the validity of: ",database,", ",input,", ",output,", and blast_params: ",blast_params," ."))}
                )
        }
        
        if(!write.only){
                
                       # default parameters that are not returned in the hit table
                       # http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
                       non_returned_params <- vector(mode = "character", length = 24)
                       non_returned_params <- c("db","query_loc","outfm","out","subject","subject_loc",
                                                "show_gis","num_descriptions","num_alignments","max_target_seqs",
                                                "html","gilist","negative_gilist","entrez_query","culling_limit",
                                                "best_hit_overhang","best_hit_score_edge","dbsize","searchsp",
                                                "import_search_strategy","export_search_strategy","parse_deflines",
                                                "num_threads","remote")
        
                       split_params <- unlist(strsplit(blast_params," "))
                       param_names <- split_params[sapply(split_params,grepl,pattern="-[a-zA-Z]",perl = TRUE)]  
                       param_names <- stringr::str_replace(param_names, pattern = "-", replacement = "")
        
                       default_params <- vector(mode = "character")
                       default_params <- c("query_id","subject_id","subject_taxonomy","subject_kingdom", "perc_identity",
                                           "num_ident_matches","alig_length","mismatches", "gap_openings","q_start",
                                           "q_end","s_start","s_end","evalue","bit_score","score_raw",
                                           "query_coverage_per_subj")
        
        
                       if(taxonomy == FALSE)
                               default_params <- default_params[!((default_params == "subject_taxonomy") | (default_params == "subject_kingdom"))]
        
        
                       additional_params <- param_names[(!is.element(param_names,default_params)) & (!is.element(param_names,non_returned_params))]
        
                       colNames <- c(default_params,additional_params)
       
                       # define the colClasses for faster file streaming
                       if(taxonomy == TRUE)
                               col_Classes <- c(rep("character",4), rep("numeric",13))
        
                       if(taxonomy == FALSE)
                               col_Classes <- c(rep("character",2), rep("numeric",13))
        
                       
                       tryCatch(
                                {
                                  if(!sql_database){
                                
                                               hit_table <- data.table::fread(input = output, sep = "\t", 
                                                                             header = FALSE, colClasses = col_Classes)
 
                                               data.table::setnames(hit_table, old = paste0("V",1:length(colNames)),
                                                                    new = colNames)
                        
                                               data.table::setkeyv(hit_table, c("query_id","subject_id"))
                                    
                                               return(hit_table) 
                        
                                   }
                        
                                   if(sql_database){
                                
                                        
                                               blast_sql_db <- dplyr::src_sqlite(paste0("_blast_db",f_sep,"blast_sql_db.sqlite3"), create = TRUE)
                                
                                               connect_db <- RSQLite::dbConnect("SQLite", dbname = paste0("_blast_db",f_sep,"blast_sql_db.sqlite3"))
                                
                                               RSQLite::dbWriteTable(connect_db, name = "hit_tbl",value = output,row.names = FALSE,
                                                                      header = FALSE, sep = "\t", overwrite = TRUE)
                                
                                               blast_sqlite <- dplyr::tbl(dplyr::src_sqlite(paste0("_blast_db",f_sep,"blast_sql_db.sqlite3")),"hit_tbl")
                                
                                               # dplyr::rename(blast_sqlite,colNames)
                                
                                               return(blast_sqlite)
                                
                                    }
                          
                        
                               }, error = function(){ stop(paste0("File ",output," could not be read properly, please check whether BLAST ",
                                  "correctly wrote a resulting BLAST hit table to ",output," ."))}
                       )
        }
}



#' @title Advanced interface function to makeblastdb
#' @description This function provides a simple, but powerful interface
#' between the R language and 'makeblastdb'.
#' @param database_file a character string specifying the path to the input file
#' that shall be transformed to a blast-able database.
#' @param params a character string specifying the arguments in the same notation as
#' calling makeblastdb from a shell like environment that shall be handed to the 'makeblastdb' call.
#' Examples could be: \code{params} = "-input_type fasta -dbtype prot -hash_index".
#' @param folder a character string specifying the internal folder in which the database shall be stored. Default is \code{folder} = "_blast_db".
#' @param path a character string specifying the path to the makeblastdb program (in case you don't use the default path). 
#' Default is \code{path} = \code{NULL}.
#' @references
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. & Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., & Schäffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' 
#' # make the A. thaliana genome to a blast-able database
#' advanced_makedb( database_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'                  params = "-input_type fasta -dbtype prot -hash_index" )
#'                  
#'                  
#' }
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_blast}}, \code{\link{advanced_blast}}
#' @export
advanced_makedb <- function(database_file, params, folder = "_blast_db/", path = NULL){
        
        
        outfile <- set_path(database_file, add.folder = folder)
        
        if(!file.exists(folder))
                dir.create(folder)
        
        file.copy(database_file,folder)
        
        f_sep <- .Platform$file.sep
        input_name <- unlist(strsplit(outfile,f_sep))
        input_name <- input_name[length(input_name)]

        tryCatch({
                
                if(is.null(path))
                        system(paste0("makeblastdb -in ",paste0(folder,input_name)," ",params))
                
                if(!is.null(path))
                        system(paste0("export PATH=makeblastdb -in ",paste0(folder,input_name)," ",params))
                
        }, error = function(){ stop(paste0("Something went wrong with the makeblastdb call.","\n",
                               "Please check your aruments: ",params," and database_file: ",database_file))}
                
                )
}





