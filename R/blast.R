#' @title Interface function to BLAST+
#' @description This function performs a BLAST+ search of a given set of sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @references 
#' 
#' http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
#' 
#' http://blast.ncbi.nlm.nih.gov/Blast.cgi
#' @examples \dontrun{
#' # performing a BLAST search using blastp (default)
#' blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'       subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' 
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
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_best}}
#' @export
blast <- function(query_file, subject_file, format = "fasta",
                  blast_algorithm = "blastp", eval = "1E-5",
                  path = NULL, comp_cores = 1, blast_params = NULL){
        
        if(!is.element(blast_algorithm,c("blastp")))
                stop("Please choose a valid BLAST mode.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # initialize the BLAST search
        query.dt <- set_blast(file = query_file)[[1]]
        # make a BLASTable databse of the subject
        database <- set_blast(file = subject_file, makedb = TRUE)[[2]]
        # create an internal folder structure for the BLAST process 
        input = paste0("_blast",f_sep,"blastinput.fasta") 
        output = paste0("_blast",f_sep,"blastresult.csv")
        
        
        if(!file.exists(paste0("_blast",f_sep))){
                
                dir.create("_blast")
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


#' @title Function to perform a BLAST best hit search
#' @description This function performs a blast+ search (best hit) of a given set of protein sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". Default is "fasta".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn",... .
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details Given a set of protein sequences A, a best hit blast search is being performed from A to database.
#' @examples \dontrun{
#' 
#' # performing gene orthology inference using the best hit (BH) method
#' blast_best(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
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
#' @seealso \code{\link{blast}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{set_best}}
#' @export
blast_best <- function(query_file, subject_file, format = "fasta", 
                       blast_algorithm = "blastp", eval = "1E-5",
                       path = NULL, comp_cores = 1, blast_params = NULL){
        
        # default parameters for best hit filtering
        default_pars <- "-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1"
        
        
        # performing a BLAST search from query against subject: blast(query,subject)
        # using the BLAST parameter: '-max_target_seqs 1' allows to retain
        # only the best hit result and therefore, speeds up the BLAST search process
        hit_tbl.dt <- blast(query_file = query_file, 
                            subject_file = subject_file, 
                            format = format, path = path, 
                            comp_cores = comp_cores,blast_params = ifelse(!is.null(blast_params),
                                                                          paste0(blast_params,
                                                                          " ",default_pars),
                                                                          default_pars))
     
        #select only the best hit (E-value) 
        #besthit_tbl <- hit_tbl.dt[ , list(.SD[ , subject_id], lapply(.SD[ , evalue],min)),
        #                         by = key(hit_tbl.dt)]        
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


#' @title Function to perform a BLAST reciprocal best hit (RBH) search
#' @description This function performs a blast+ search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
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
#'
#' @examples \dontrun{
#' # performing gene orthology inference using the best hit (BH) method
#' blast_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'           subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
#' 
#' # use multicore processing
#' blast_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
#'            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'            comp_cores = 2)
#' }
#'
#' @return A data.table as returned by the blast() function, storing the geneids
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @seealso \code{\link{blast}}, \code{\link{blast_best}}, \code{\link{advanced_blast}}, \code{\link{set_best}}
#' @export
blast_rec <- function(query_file, subject_file, format = "fasta", 
                      blast_algorithm = "blastp", eval = "1E-5",
                      path = NULL, comp_cores = 1, blast_params = NULL){
        
        orthoA <- blast_best(query_file,subject_file, 
                             format = format, blast_algorithm = blast_algorithm,
                             path = path, comp_cores = comp_cores, blast_params = blast_params)
        
        orthoB <- blast_best(subject_file,query_file, 
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
#' @description This function reads a cds fasta file using read.cds(), translates each
#' sequence into the corresponding aminoacid sequence and is able to create a blast
#' database from that.
#' @param file a character string specifying the path to the file storing the cd.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param makedb TRUE or FALSE whether a database should be created or not.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @return A list with [[1]] a data.table storing the gene id in the first column and
#' the corresponding dna and aminoacid sequence as string in the second and third column
#' respectively and [[2]] the name of the database created.
#' @examples \dontrun{
#'  # running the set function to see an example output
#'  head(set_blast(file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'))[[1]] , 2)
#' }
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{blast}}
#' @export
set_blast <- function(file, format = "fasta", makedb = FALSE, path = NULL, ...){
        
        # HERE WE NEED TO INSERT SOME QUALITY CONTROL
        # CONCERNING THE CASE THAT A POTENTIAL USER
        # COULD INSERT SOME NON-CDS SEQUENCES
        
        # read cds file
        # copy the data.table due to this discussion:
        # http://stackoverflow.com/questions/8030452/pass-by-reference-the-operator-in-the-data-table-package
        dt <- data.table::copy(read.cds(file = file, format = "fasta", ...))
        
        if(!is.data.table(dt))
                stop("Your CDS file was not corretly transformed into a data.table object.")
        
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        
        # WE COULD ALSO INSERT A IF-STATEMENT CHECKING THAT
        # IN CASE THE USER INSERTS AN AA-FASTA FILE,
        # THE TRANSLATION PROCESS IS BEING OMITTED
        
        
        # When using data.tables within packaes, always make sure
        # 'data.table' is included in the DESCRIPTION file as 'Imports' AND
        # in the NAMESPACE file with 'imports(data.table)' -> @import when using roxygen2
        # http://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package
        # https://github.com/hadley/dplyr/issues/548
        
        # translate dna to aa sequences
        #dt[ , aa := as.vector(sapply(seqs,transl)), by = geneids]
        
        tryCatch(
          {
                     dt[ , aa := transl(seqs), by = geneids]
                     
          }, error = function() {stop(paste0("The input coding sequences could not be translated properly to amino acid sequences.",
          "\n Please check whether ",file, " stores valid coding sequences."))}
        )

        # makedb
        dbname <- vector(mode = "character", length = 1)
        filename <- unlist(strsplit(file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
        filename <- filename[length(filename)]
        
        
        if(makedb){
                if(!file.exists(paste0("_database",f_sep))){ 
                        
                        dir.create("_database") 
                
                }
                
                dbname <- paste0("_database",f_sep,"out_",filename,"_translate.fasta")
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
                                               " -input_type fasta -dbtype prot")
                                       )
                        } else {
                                system(
                                        paste0("export PATH=$PATH:",path,"; makeblastdb -in ",
                                               dbname," -input_type fasta -dbtype prot")
                                      )
                        }
                }, error = function(){ stop(paste0("makeblastdb did not work properly. The default parameters are: \n",
                "-input_type fasta -dbtype prot . \n","Please check that you really want to work with a protein database.\n",
                "Additionally check: ",dbname," ."))}
                )
        }

        return(list(dt, dbname))
}



#' @title Advanced interface function to BLAST+
#' @description This function performs a BLAST+ search of a given set of sequences and a given set of parameters against a BLAST database.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk".
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, e.g. "blastp","blastn","tblastn","deltablast",... .
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param taxonomy a logical value specifying whether the subject taxonomy shall be returned based on NCBI Taxonomy queries. Default is \code{FALSE}.
#' @param taxdb_path a character string specidying the path to the lical taxonomy database.
#' @details
#'  
#' Following BLAST programs and algorithms can be assigned to \code{blast_algorithm}:
#' 
#' "blastp", "blastn", "megablast","psi-blast", "phi-blast", "delta-blast", "blastx", "tblastn", "tblastx"
#' 
#' @author Hajk-Georg Drost and Sarah Scharfenberg
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
#'       
#'              
#' }
#'
#' @return the \code{advanced_blast} function creates a folder named '_blast' and stores
#' the hit table returned by BLAST in this folder.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_best}}
#' @export
advanced_blast <- function(query_file, subject_file, format = "fasta", 
                           blast_algorithm = "blastp", path = NULL,
                           blast_params = NULL,taxonomy = FALSE,
                           taxdb_path = NULL){
        
        # http://blast.ncbi.nlm.nih.gov/Blast.cgi
        if(!is.element(blast_algorithm,c("blastp","blastn", "megablast",
                                         "psiblast", "phiblast", "deltablast",
                                         "blastx","tblastn","tblastx")))
                stop("Please choose a valid BLAST mode.")
        
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # initialize the BLAST search
        
        query.dt <- set_blast(file = query_file, format = format)[[1]]
        # make a BLASTable databse of the subject
        database <- set_blast(file = subject_file, format = format, makedb = TRUE)[[2]]
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # create an internal folder structure for the BLAST process 
        input = paste0("_blast",f_sep,"blastinput.fasta") 
        output = paste0("_blast",f_sep,"blastresult.csv")
        
        
        if(!file.exists(paste0("_blast",f_sep))){
                
                dir.create("_blast")
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
                                
                                               if(!is.null(taxdb_path)){
                                                       system( paste0("export BLASTDB=$BLASTDB:",taxdb_path,"; ",blast_algorithm," -db ",database," -query ",input,
                                                               " -out ", output ," ",blast_params,
                                                               " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                               length mismatch gapopen qstart qend sstart send evalue 
                                                               bitscore score qcovs'")
                                                              )
                                                }
                                        }
                              }, error = function(){ stop(paste0("taxdb could not be included to the BLAST search. \n",
                                                     "Please check the validity of the path: ",taxdb_path," .\n",
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
                        
                                         if(!is.null(taxdb_path)){     
                                
                                                 system(
                                                        paste0("export PATH=$PATH:",path,"; ","export BLASTDB=$BLASTDB:",taxdb_path,"; ",blast_algorithm," -db ",
                                                        database," -query ",input," -out ", output ," ",blast_params,
                                                        " -outfmt '6 qseqid sseqid staxids sskingdoms pident nident 
                                                        length mismatch gapopen qstart qend sstart send evalue bitscore 
                                                        score qcovs'")
                                                       )
                                         }
                       
                                  }
                         }, error = function(){stop(paste0("taxdb could not be included to the BLAST search. \n",
                                                            "Please check the validity of the path: ",taxdb_path," .\n",
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
                        hit_table <- data.table::fread(input = output, sep = "\t", 
                                                       header = FALSE, colClasses = col_Classes)
 
                        data.table::setnames(hit_table, old = paste0("V",1:length(colNames)),
                                             new = colNames)
        
                        data.table::setkeyv(hit_table, c("query_id","subject_id"))
        
                        return(hit_table)      
                        
                }, error = function(){ stop(paste0("File ",output," could not be read properly, please check whether BLAST ",
                "correctly wrote a resulting BLAST hit table to ",output," ."))}
        )
        
}
