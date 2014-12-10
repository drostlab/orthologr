#' @title Function to perform a global alignment search for orthologous genes.
#' @description This function performs a protein alignment for each sequences
#' from the query_file to each of sequence of the subject_file and returns
#' a data.table containing the ids of the best global alignment pair and the corresponding
#' score parameter.
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_file a character string specifying the path to the sequence file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein". In case of "cds", sequence are translated to protein sequences,
#' Default is \code{seq_type} = "protein".
#' @param format a character string specifying the file format of the sequence file. Default is "fasta".
#' @param tool a character string specifying the global alignment tool that shall be used, e.g. "ggsearch","ssearch".
#' @param path a character string specifying the path to the tool executable (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores to be used for.
#' @param statistics a boolean value specifying whether the statistics should be printed and appear in the hit table.
#' @param details a boolean value specifying whether a detailed hit table shall be returned or only the evalue of the corresponding
#' ortholog pairs.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param quiet a boolean value suppressing any output if true
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @examples \dontrun{
#' 
#' # performing a search between two protein fastas using Needleman Wunsch algorithm
#' alignmentSearch( query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'                  subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'))
#'
#' # performing a search between two cds fastas in multicore mode
#' alignmentSearch( query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                  subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'), 
#'                  comp_cores=2, seq_type = "cds")
#'                  
#'                  
#' # performing a search between two cds fastas using Smith Waterman algorithm
#' alignmentSearch( query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#'                  subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'), 
#'                  tool="ssearch",comp_cores=1, seq_type = "cds")
#'  
#' # performing a search between two protein fastas getting the most detailed output
#' alignmentSearch( query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'                  subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'                  path="/home/sarah/Programs/fasta-36.3.7/bin/", details=TRUE,
#'                  statistics=TRUE )
#'    }
#'@return a data.table containing query id and subject id as well as a series of score values.
#'@export 
alignmentSearch <- function(query_file, subject_file, seq_type = "prot",
                            format = "fasta", tool = "ggsearch",
                            path = NULL, comp_cores = 1, details=FALSE,
                            statistics = FALSE, clean_folders = FALSE, quiet=TRUE){      
        
        if(!is.alignment_search_tool(tool)){
                stop("Please choose a global alignment tool that is supported by this function.")
        }
        
        
        f_sep <- .Platform$file.sep
        options <- ""
        
        if(!file.exists(paste0("_global_orthologs",f_sep))){ 
                
                dir.create("_global_orthologs") 
                
        }
                
        if(comp_cores > parallel::detectCores())
                stop("You assigned more cores to the comp_cores argument than are availible on your machine.")
        multicore <- (comp_cores > 1)
        
        if(seq_type == "cds"){
                
                # read cds
                q_dt <- data.table::copy(read.cds(file = query_file, format = "fasta"))   
                # omit empty sequences
                q_dt <- q_dt[,.SD[sapply(seqs,function(x){return(! (is.na(x) || x=="") )})]]
                # omit sequences taht are not multiples of 3
                q_dt <- q_dt[,.SD[sapply(seqs,function(x){return(nchar(x)%%3==0)})]]
                # omit sequences consisting of others than ACGT
                q_dt <- q_dt[,.SD[sapply(seqs,is.dnaSequence)]]
                # translate query
                tryCatch(
                { q_dt[ , aa := transl(seqs), by = geneids]
                 }, error = function() {stop(paste0("The input coding sequences could not be translated properly to amino acid sequences.",
                                           ,"\n"," Please check whether ",file, " stores valid coding sequences."))}
                )
                
                # write protein sequences to library
                filename <- unlist(strsplit(query_file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                queries <- paste0("_global_orthologs",f_sep,"out_",filename,"_translate.fasta")
                seqinr::write.fasta(as.list(q_dt[ , aa]),
                                    names = q_dt[ , geneids],
                                    nbchar = 80, open = "w",
                                    file.out = queries 
                )
                
                n_seq <- nrow(q_dt)
                
                s_dt <- read.cds(file = subject_file, format = "fasta")
                # omit empty sequences
                s_dt <- s_dt[,.SD[sapply(seqs,function(x){return(! (is.na(x) || x=="") )})]]
                # omit sequences taht are not multiples of 3
                s_dt <- s_dt[,.SD[sapply(seqs,function(x){return(nchar(x)%%3==0)})]]
                # omit sequences consisting of others than ACGT
                s_dt <- s_dt[,.SD[sapply(seqs,is.dnaSequence)]]
                # translate subject
                tryCatch(
                { s_dt[ , aa := transl(seqs), by = geneids]
                }, error = function() {stop(paste0("The input coding sequences could not be translated properly to amino acid sequences.",
                                                   ,"\n"," Please check whether ",file, " stores valid coding sequences."))}
                                )
                
                # write protein sequences to library
                filename <- unlist(strsplit(subject_file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                dbname <- paste0("_global_orthologs",f_sep,"out_",filename,"_translate.fasta")
                seqinr::write.fasta(as.list(s_dt[ , aa]),
                                    names = s_dt[ , geneids],
                                    nbchar = 80, open = "w",
                                    file.out = dbname 
                )
                
        }
             
        if(seq_type == "protein"){
                n_seq <- length(Biostrings::readAAStringSet(filepath = query_file, format = "fasta"))
                queries <- query_file
                dbname <- subject_file
        }

        # search
        options <- paste0(options,"-p -T ",comp_cores)
         
        if( tool == "ggsearch"){
                ggsearch(file=queries, library_file=dbname, path = path, options = options,
                         parse_output_to = "_global_orthologs/search.out")
        }
        if( tool == "ssearch"){
                ssearch(file=queries, library_file=dbname, path = path, options = options,
                         parse_output_to = "_global_orthologs/search.out")
        }
        # parse output file
        con  <- file("_global_orthologs/search.out", open = "r")
        collect <- list()
        new_names <- c("query_id","subject_id","e_value")
        current.line <- 1
        while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
                
                # grep >>> begin new entry
                if(grepl(">>>",line)){
                        
                        line <- gsub("^\\s+","",line)
                        splits <- strsplit(line, ">>>| ")
                        query_id <- splits[[1]][2]
                        end=FALSE
                        while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0 && !end){
                        
                                if(statistics && grepl("^Statistics:",line)){
                                        splits <- strsplit(line, " |=")
                                        mu <- splits[[1]][9]
                                        var <- splits[[1]][12]
                                        if(!quiet){
                                                print(query_id)
                                                for(i in c(1:5)){
                                                        print(line)
                                                        line <- readLines(con, n = 1, warn = FALSE)
                                                }                                
                                        }
                                }
                                
                                if(grepl(">>",line)){
#                                 if(grepl("The best scores are:",line)){
#                                         line <- readLines(con, n = 1, warn = FALSE)
#                                         splits <- strsplit(line, " ")
#                                         subject_id <- splits[[1]][1]
#                                         e_val <- splits[[1]][length(splits[[1]])]
                                        splits <- strsplit(line,">>| ")
                                        subject_id <- splits[[1]][2]
                                        line <- readLines(con, n = 1, warn = FALSE)
                                        splits <- strsplit(line,":")
                                        z_score <- strsplit(splits[[1]][3],"( )+")[[1]][2]
                                        bit_score <- strsplit(splits[[1]][4],"( )+")[[1]][2]
                                        e_val <- strsplit(splits[[1]][5],"( )+")[[1]][2]
                                        
                                        if(details){
                                                if(statistics){
                                                        collect[[current.line]] <- c(query_id,subject_id, e_val, bit_score, z_score, mu, var)
                                                        new_names <- c("query_id","subject_id","e_value","bit_score","Z_Score","mu","var")
                                                }
                                                else{
                                                        collect[[current.line]] <- c(query_id,subject_id, e_val, bit_score, z_score)
                                                        new_names <- c("query_id","subject_id","e_value","bit_score","Z_Score")
                                                }
                                        }
                                        else{
                                                if(statistics){
                                                        collect[[current.line]] <- c(query_id,subject_id, e_val, mu, var)
                                                        new_names <- c("query_id","subject_id","e_value","mu","var")
                                                }
                                                else{
                                                        collect[[current.line]] <- c(query_id,subject_id, e_val)
                                                        new_names <- c("query_id","subject_id","e_value")
                                                }
                                                
                                        }
                                        
                                        current.line <- current.line + 1
                                        query_id <- NULL
                                        subject_id <- NULL
                                        end=TRUE
                                }
                        }
                }
        } 
        close(con)
        
        hit.table <- data.table(do.call(rbind,collect))
        old_names <- names(hit.table)
        setnames(hit.table, old=old_names, new = new_names )
        
        if(clean_folders){
                clean_all_folders("_global_orthologs")
        }
        
        return(hit.table)
}
