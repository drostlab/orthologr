#' @title Perform an advanced BLAST+ search
#' @description This function performs a BLAST+ search of a given set of sequences and a given set of parameters against a BLAST database.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' In case a NCBI database is used, it is possible to specify \code{subject_file} = \code{"nr"}, \code{"plaza"}, \code{"cdd_delta"}. Other databases
#' will be added in future versions of \pkg{orthologr}. When specifying \code{subject_file} = "nr", "plaza", ...
#' make sure the corresponding database ("nr", "plaza","cdd_delta" ...) is stored in the blast working directory \code{"_blast_db"}.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: \code{seq_type} = \code{"cds"}, \code{seq_type} = \code{"protein"}, or \code{seq_type} = \code{"dna"}.
#' In case of \code{seq_type} = \code{"cds"}, sequence are translated to protein sequences,
#' in case of \code{seq_type} = \code{"dna"}, cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = \code{"cds"}.
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format} = \code{"fasta"}, \code{format} = \code{"gbk"}.
#' @param blast_algorithm a character string specifying the BLAST algorithm that shall be used, 
#' e.g. \code{blast_algorithm} = \code{"blastp"}, \code{blast_algorithm} = \code{"blastn"}, \code{blast_algorithm} = \code{"tblastn"},
#' \code{blast_algorithm} = \code{"deltablast"}.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param blast_params a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running BLAST.
#' @param makedb_type a character string specifying the sequence type stored in the BLAST database
#' that is generated using 'makeblastdb'. Options are: "protein" and "nucleotide". Default is \code{makedb_type} = \code{"protein"}.
#' @param taxonomy a logical value specifying whether the subject taxonomy shall be returned based on NCBI Taxonomy queries. 
#' Default is \code{taxonomy} = \code{FALSE}.
#' @param db_path a character string specidying the path to the local BLAST database.
#' @param sql_database a logical value specifying whether an SQL database shall be created to store the BLAST output.
#' This is only useful when the amount of data does not fit in-memory anymore.
#' @param session_name a character string specifying the name of the BLAST output. Default is \code{session_name} = \code{NULL}.
#' @param write.only a logical value specifying whether the resulting blast output should be returned as
#' data.table or database connection. Default is \code{write.only} = \code{FALSE}.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @details
#' 
#' The \code{advanced_blast} function allows you to use the BLAST+ stand-alone specific notation
#' within the R environment. The \code{blast_params} argument allows the user to specify BLAST+ specific command line
#' parameters and the resulting hit table is then returned as \pkg{data.table} object.
#' 
#' Optionally the \code{sql_database = TRUE} argument allows you to store the hit table in an SQLite database
#' using the flexible query notation provided by the \pkg{dplyr} (>= 0.3) package.
#' 
#' The \code{advanced_blast} function works very similar to the command line version of BLAST+,
#' so when specifying a \code{blast_algorithm} such as \emph{blastp}, \emph{tblastn}, \emph{deltablast}, etc.
#' the corresponding databases must be stored and accessible accoring to the BLAST stand-alone notation. 
#' 
#' 
#' Following BLAST programs and algorithms can be assigned to \code{blast_algorithm}:
#' 
#' \itemize{
#' 
#' \item "blastp"
#' \item "blastn"
#' \item "megablast"
#' \item "psiblast"
#' \item "phiblast"
#' \item "deltablast"
#' \item "blastx"
#' \item "tblastn"
#' \item "tblastx"
#' 
#' }
#' 
#' The intention of this function is to provide the user with an easy to use interface function
#' to BLAST+.
#' 
#' When using \code{taxonomy} =  \code{TRUE}, then the taxdb.btd and taxdb.bti must be stored
#' in the same folder as BLAST-able database. In case you need to specify the path to your BLAST database,
#' use the \code{db_path} argument.
#' 
#' The \code{advanced_blast} function can also store the BLAST output CSV file in a SQLite database.
#' This works only with dplyr version >= 0.3 . To store the BLAST output in an SQLite database and
#' to receive an SQLite connection as specified by \code{tbl} in \code{dplyr} please use the \code{sql_database} = \code{TRUE} argument. 
#' 
#' In case you want to use deltablast, you need to all cdd_delta.* files into the folder "_blast_db".
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
#' BLAST Command Line Applications User Manual: \url{http://www.ncbi.nlm.nih.gov/books/NBK1763/}
#' 
#' Blast Program Selection Guide: \url{http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide}
#' 
#' BLAST web interface: \url{http://blast.ncbi.nlm.nih.gov/Blast.cgi}
#' 
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. and Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' Gish, W. and States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L., and Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W., and Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., and Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#'
#' Zhang, J. and Madden, T.L. (1997) "PowerBLAST: A new network BLAST application for interactive or automated sequence analysis and annotation." Genome Res. 7:649-656.
#'
#' Morgulis A., Coulouris G., Raytselis Y., Madden T.L., Agarwala R., and Schaeffer A.A. (2008) "Database indexing for production MegaBLAST searches." Bioinformatics 15:1757-1764.
#'
#' Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., and Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
#' 
#' 
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
#'       seq_type = "protein", blast_algorithm = "blastp", 
#'       blast_params = "-evalue 1E-5 -num_threads 2")       
#'  
#' 
#'              
#' # DELTA-BLAST
#' #
#' # you can also use deltablast to perform BLAST searches
#' # make sure you have the cdd_deltablast.* files stored in the folder "_blast_db"                          
#' advanced_blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                blast_algorithm = "deltablast", blast_params = "-evalue 1E-5 -num_threads 2")                                       
#'                              
#'                                                                  
#'                                                                                                                                          
#' # when performing an advanced BLAST search, you can easily select the best hit using
#' library(dplyr)
#' 
#' advB <- advanced_blast(
#'            query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'            subject_file = system.file('seqs/ortho_lyra_cds_1000.fasta', package = 'orthologr'),
#'            seq_type = "cds",blast_algorithm = "blastp", 
#'            blast_params = "-evalue 1E-5 -num_threads 1")
#'                                  
#' best_hit <- advB %>% group_by(query_id) %>% summarise(min(evalue))
#'
#'
#' # using a SQLite database to store the BLAST output and
#' # select some example data
#' 
#' library(dplyr)
#' 
#' sqlE <- advanced_blast(
#'              query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'              subject_file = system.file('seqs/ortho_lyra_cds_1000.fasta', package = 'orthologr'),
#'              seq_type = "cds",blast_algorithm = "blastp", 
#'              blast_params = "-evalue 1E-5 -num_threads 1", sql_database = TRUE)
#'
#' head(sqlE)
#' 
#' glimpse(sqlE)
#' 
#' # select all rows that have an evalue of zero
#' filter(sqlE,evalue == 0)
#' 
#' # select the best hit using the evalue criterion
#' sqlE %>% group_by(query_id) %>% summarise(best_hit_eval = min(evalue))
#' 
#'                                                  
#'                                                                                                   
#' # using advanced_blast() with external BLAST path                                                                                                                                                   
#' advanced_blast(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#'                blast_algorithm = "blastp", blast_params = "-evalue 1E-5 -num_threads 2",
#'                path = "path/to/blastp/")
#' 
#' }
#'
#' @return the \code{advanced_blast} function creates a folder named '_blast_db' and stores
#' the hit table returned by BLAST in this folder.
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{blast}}, \code{\link{set_blast}}, \code{\link{advanced_makedb}}
#' @export
advanced_blast <- function(query_file, 
                           subject_file, 
                           seq_type        = "cds",
                           format          = "fasta", 
                           blast_algorithm = "blastp", 
                           path            = NULL,
                           blast_params    = NULL, 
                           makedb_type     = "protein",
                           taxonomy        = FALSE, 
                           db_path         = NULL,
                           sql_database    = FALSE, 
                           session_name    = NULL,
                           write.only      = FALSE, 
                           clean_folders   = FALSE){
        
        # http://blast.ncbi.nlm.nih.gov/Blast.cgi
        if(!is.element(blast_algorithm,c("blastp","blastn", "megablast",
                                         "psiblast", "phiblast", "deltablast",
                                         "blastx","tblastn","tblastx")))
                stop("Please choose a valid BLAST mode.")
        
        if(is.element(subject_file,c("nr","plaza","cdd_delta"))){
                
                use_ncbi_database <- TRUE
                
        } else {
                
                use_ncbi_database <- FALSE
                
        }
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        aa <- geneids <- NULL
        
        # initialize the BLAST search
        
        query.dt <- set_blast( file     = query_file, 
                               format   = format, 
                               seq_type = seq_type )[[1]]
        
        if(!use_ncbi_database){
                
                # make a BLASTable databse of the subject
                database <- set_blast( file        = subject_file, 
                                       format      = format, 
                                       seq_type    = seq_type, 
                                       makedb      = TRUE, 
                                       makedb_type = makedb_type )[[2]]
        }
        
        if(use_ncbi_database){
                
                if(use_ncbi_database == "nr")
                        database <- "nr" 
                
                if(use_ncbi_database == "plaza")
                        database <- "plaza"
                
                if(use_ncbi_database == "cdd_delta")
                        database <- "cdd_delta" 
                
        }
        
        if(is.null(session_name))
                session_name <- get_filename(query_file)
        
        # create an internal folder structure for the BLAST process 
        input = paste0(session_name,"_blastinput.fasta") 
        output = paste0(session_name,"_blastresult.csv")
        
        
        if(!file.exists(file.path(tempdir(),"_blast_db"))){
                
                dir.create(file.path(tempdir(),"_blast_db"))
        }
        
        currwd <- getwd()
        setwd(file.path(tempdir(),"_blast_db"))
        
        # write query fasta file 
#         write_AA <- as.list(query.dt[ ,aa])
#         name <- query.dt[ , geneids]      
        
        tryCatch({
                
                seqinr::write.fasta( sequences = as.list(query.dt[ ,aa]), 
                                     names     = query.dt[ , geneids],
                                     nbchar    = 80,
                                     open      = "w",
                                     file.out  = input )
        
        }, error = function(e){ stop("File ",input," could not be written properly to the internal folder environment.","\n",
                                   "Please check: ",query_file, " and ",subject_file,".")}
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
                }, error = function(e){ stop(paste0("taxdb could not be included to the BLAST search. \n",
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

        }, error = function(e){stop("The advanced BLAST search did not work properly.","\n",
                                  "Please check the validity of: ",database,", ",input,", ",output,", and blast_params: ",blast_params," .")}
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
        
        }, error = function(e){stop("taxdb could not be included to the BLAST search. \n",
                                          "Please check the validity of the path: ",db_path," .\n",
                                          "Additionally, check the validity of: ",database,", ",input,",
                                          ",output,", and blast_params: ",blast_params," .")}
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
        
        }, error = function(e){stop("The advanced BLAST search did not work properly.","\n",
                                  "Please check the validity of: ",database,", ",input,", ",output,", and blast_params: ",blast_params," .")}
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
                
                hit_table <- data.table::fread(input      = output, 
                                               sep        = "\t", 
                                               header     = FALSE,
                                               colClasses = col_Classes)
                
                data.table::setnames(hit_table, old = paste0("V",1:length(colNames)),
                                     new = colNames)
                
                data.table::setkeyv(hit_table, c("query_id","subject_id"))
                
                # return to the global working directory
                setwd(file.path(currwd))
                
                if(clean_folders)
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                
                return(hit_table) 
                
        }
        
        if(sql_database){
                
                
                blast_sql_db <- dplyr::src_sqlite("blast_sql_db.sqlite3", create = TRUE)
                
                connect_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = "blast_sql_db.sqlite3")
                
                
                DBI::dbWriteTable(connect_db, 
                                  name      = "hit_tbl",
                                  value     = output,
                                  row.names = FALSE,
                                  header    = FALSE, 
                                  sep       = "\t", 
                                  overwrite = TRUE)
                
                on.exit({
                        
                        DBI::dbDisconnect(connect_db)
                        
                })
                
                blast_sqlite <- dplyr::tbl(dplyr::src_sqlite("blast_sql_db.sqlite3"),"hit_tbl")
                
#               dplyr::rename(blast_sqlite,)

                # return to the global working directory
                setwd(file.path(currwd))
                
                if(clean_folders){
                        
                        warnings("Are you sure you want to clean all folders? The SQLite database
                                 you want to connect with will be deleted as well.")
                        clean_all_folders(file.path(tempdir(),"_blast_db"))
                        
                        cat("\n")
                        warnings("Database has been deleted.")
                }
                
                return(blast_sqlite)
                
        }
        
        
        }, error = function(e){ stop("File ",output," could not be read properly, please check whether BLAST ",
                                   "correctly wrote a resulting BLAST hit table to ",output," .")}
        )
} 

if(write.only){
        
        stop("Not implemented yet.")
        
}


}






