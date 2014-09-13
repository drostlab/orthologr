#' @title Interface function to BLAST+
#' @description This function performs a blast+ search of a given set of protein sequences against a given database.
#' @param query.dt a data.table returned as first list element of the set() function.
#' @param database a string specifying the path to the BLAST database, as returned
#' as second list element of the set() function.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param input a character string specifying the path to the input file.
#' @param output a character string specifying the path to the output file.
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run BLAST searches
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @examples \dontrun{
#' # performing a best hit BLAST search
#' blast( set("data/ortho_thal_cds.fasta")[[1]] , set("data/ortho_lyra_cds.fasta", makedb = TRUE)[[2]])
#' 
#' # in case you are working with a multicore machine, you can also run parallel
#' # BLAST computations using the comp_cores parameter: here with 2 cores
#' blast( set("data/ortho_thal_cds.fasta")[[1]] , 
#'        set("data/ortho_lyra_cds.fasta", makedb = TRUE)[[2]],
#'        comp_cores = 2)
#' }
#'
#' @return A data.table storing the BLAST hit table returned by the BLAST program.
#' @export
blast <- function(query.dt, database, eval = "1E-5",
                  input = "_blast/_blastinput.fasta", 
                  output = "_blast/_blastresult.csv",
                  path = NULL, comp_cores = 1){
        
        if(!file.exists("_blast/")){
                
                dir.create("_blast")
        }
        
        # here to come: a parallelized version of the BLAST process
        
        # determine the number of cores on a multicore machine
        cores <- parallel::makeForkCluster(parallel::detectCores(all.tests = FALSE, logical = FALSE))
        
        # in case one tries to use more cores than are available
        if(comp_cores > cores)
                stop("You chose more cores than are available on your machine.")
        
        
        # DO NOT RUN THIS with big data, as it tries to blast all sequences at a time
        write_AA <- as.list(query.dt[ ,aa])
        name <- query.dt[ ,geneids]      
                
        seqinr::write.fasta(write_AA, names = name,
                            nbchar = 80,open = "w",
                            file.out = input)
                
        if(is.null(path)){
                system(
                        paste0("blastp -db ",database," -query ",input,
                               " -evalue ",eval," -out ", output ," -outfmt 6", 
                               " -num_threads ", comp_cores)
                )
        } else {
                system(
                        paste0("export PATH=$PATH:",path,"; blastp -db ",
                               database," -query ",input," -evalue ",
                               eval," -out ", output ," -outfmt 6", 
                               " -num_threads ", comp_cores)
                        )
        }
        
        
        # additional blast parameters can be found here:
        # http://www.ncbi.nlm.nih.gov/books/NBK1763/table/CmdLineAppsManual.T.options_common_to_al/?report=objectonly
        blast_table_names <- c("query_id","subject_id","perc_identity","alig_length","mismatches",
                               "gap_openings","q_start","q_end","s_start","s_end","evalue","bit_score")
        
        # define the colClasses for faster file streaming
        c_Classes <- c(rep("character",2), rep("numeric",10))
                
        hit_table <- data.table::fread(input = output, sep = "\t", header = FALSE, colClasses = c_Classes)
        data.table::setnames(hit_table, old = paste0("V",1:length(blast_table_names)), new = blast_table_names)
        data.table::setkey(hit_table, query_id)
        
        return(hit_table)
}


#' @title Function to perform a BLAST best hit search
#' @description This function performs a blast+ search (best hit) of a given set of protein sequences against a given database.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details Given a set of protein sequences A, a best hit blast search is being performed from A to database.
#' @examples \dontrun{
#' # performing gene orthology inference using the best hit (BH) method
#' blast_best(query_file = "data/ortho_thal_cds.fasta", subject_file = "data/ortho_lyra_cds.fasta")
#' 
#' # use multicore processing
#' blast_best(query_file = "data/ortho_thal_cds.fasta", 
#'            subject_file = "data/ortho_lyra_cds.fasta",
#'            comp_cores = 2)
#' }
#'
#' @return A data.table as returned by the blast() function, storing the geneids
#' of orthologous genes (best hit) in the first column and the amino acid sequences in the second column.
#' @export
blast_best <- function(query_file, subject_file, path = NULL, comp_cores = 1){
        
        hit_tbl.dt <- blast(query.dt = set(query_file)[[1]], 
                          database = set(subject_file, makedb = TRUE)[[2]], 
                          path = path, comp_cores = comp_cores)
        
        #uniq_geneids <- unique(hit_tbl.dt[ , query_id])
        
        besthit_tbl <- hit_tbl.dt[ , list(.SD[ , subject_id],lapply(.SD[ , evalue],min)),
                   by = key(hit_tbl.dt)]
        data.table::setkey(besthit_tbl, query_id)
        data.table::setnames(besthit_tbl, old = c("V1","V2"), new = c("subject_id","evalue"))
        
        # here we select the best hit from the data.table
        return( besthit_tbl )
              
}


#' @title Function to perform a BLAST reciprocal best hit (RBH) search
#' @description This function performs a blast+ search (reciprocal best hit) of a given set of protein sequences against a second
#' set of protein sequences and vice versa.
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for multicore BLAST computations.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details Given a set of protein sequences A and a different set of protein sequences B,
#' first a best hit blast search is being performed from A to B: blast(A,B) and afterwards
#' a best hit blast search is being performed from B to A: blast(B,A). Only protein sequences
#' that were found to be best hits in both directions are retained and returned.
#'
#'
#' @examples \dontrun{
#' # performing gene orthology inference using the best hit (BH) method
#' blast_rec(query_file = "data/ortho_lyra_cds.fasta", subject_file = "data/ortho_lyra_cds.fasta")
#' 
#' # use multicore processing
#' blast_rec(query_file = "data/ortho_thal_cds.fasta", 
#'            subject_file = "data/ortho_lyra_cds.fasta",
#'            comp_cores = 2)
#' }
#'
#' @return A data.table as returned by the blast() function, storing the geneids
#' of orthologous genes (reciprocal best hit) in the first column and the amino acid sequences in the second column.
#' @export
blast_rec <- function(query_file, subject_file, path = NULL, comp_cores = 1){
        
        orthoA <- blast_best(query_file,subject_file, path = path, comp_cores = comp_cores)
        orthoB <- blast_best(subject_file,query_file, path = path, comp_cores = comp_cores)
        return ( dplyr::semi_join(orthoA, orthoB, by = c("query_id","subject_id")) )
}


#' @title Function preparing the parameters and databases for subsequent BLAST searches
#' @description This function reads a cds fasta file using read.cds(), translates each
#' sequence into the corresponding aminoacid sequence and is able to create a blast
#' database from that.
#' @param file a character string specifying the path to the file storing the cd.
#' @param format a character string specifying the file format used to store the genome, e.g. "fasta", "gbk".
#' @param makedb TRUE or FALSE whether a database should be created or not.
#' @param ... additional arguments that are used by the seqinr::read.fasta() function.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @return A list with [[1]] a data.table storing the gene id in the first column and
#' the corresponding dna and aminoacid sequence as string in the second and third column
#' respectively and [[2]] the name of the database created.
set <- function(file, format = "fasta", makedb = FALSE, path = NULL, ...){
        
        # HERE WE NEED TO INSERT SOME QUALITY CONTROL
        # CONCERNING THE CASE THAT A POTENTIAL USER
        # COULD INSERT SOME NON-CDS SEQUENCES
        
        # read cds file
        dt <- read.cds(file=file, format="fasta", ...)
        
        # WE COULD ALSO INSERT A IF-STATEMENT CHECKING THAT
        # IN CASE THE USER INSERTS AN AA-FASTA FILE,
        # THE TRANSLATION PROCESS IS BEING OMITTED
        
        # translate dna to aa sequences (not working yet)
        dt <- dt[ , aa:=seqinr::c2s(seqinr::translate(seqinr::s2c(seqs))), by = geneids]
        # makedb
        dbname<-""
        filename <- unlist(strsplit(file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
        filename <- filename[length(filename)]
        if(makedb){
                if(!file.exists("_database/")){ 
                        
                        dir.create("_database") 
                
                }
                
                dbname <- paste0("_database/out_",filename,"_translate.fasta")
                seqinr::write.fasta(as.list(dt[ , aa]),
                                    names = dt[ , geneids],
                                    nbchar = 80, open = "w",
                                    file.out = dbname 
                                    )
                
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
        }
        return(list(dt, dbname))
}
