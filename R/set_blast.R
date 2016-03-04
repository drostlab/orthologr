#' @title Create a BLASTable database with makeblastdb
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
#'  # running the set function to see an example output
#'  head(set_blast(file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'))[[1]] , 2)
#' }
#' @seealso \code{\link{blast_best}}, \code{\link{blast_rec}}, \code{\link{advanced_blast}}, \code{\link{blast}}, \code{\link{advanced_makedb}}
#' @import data.table
#' @export

set_blast <- function(file, 
                      seq_type    = "cds",
                      format      = "fasta", 
                      makedb      = FALSE,
                      path        = NULL, 
                      makedb_type = "protein",
                      ...){
        
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
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        seqs <- aa <- geneids <- NULL
        
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
                dt <- dt[ ,.SD[sapply(seqs,function(x){return(! (is.na(x) || x=="") )})]]
                
                # omit sequences taht are not multiples of 3
                dt <- dt[ ,.SD[sapply(seqs,function(x){return(nchar(x)%%3==0)})]]
                
                # omit sequences consisting of others than ACGT
                dt <- dt[ ,.SD[sapply(seqs,is.dnaSequence)]]
                
                # translate cds to protein sequences
                tryCatch({
                        
                        dt[ , aa := transl(seqs), by = geneids]
                        
                        }, error = function(e) {stop("The input coding sequences could not be translated properly to amino acid sequences.",
                                   ,"\n"," Please check whether ",file, " stores valid coding sequences.")}
                )

        }
        
        # read file storing protein sequence -> no translation needed
        if(seq_type == "protein"){
                
                dt <- data.table::copy(read.proteome(file = file, format = "fasta", ...))
                data.table::setnames(dt, old = "seqs", new = "aa")
        }
        
        if(seq_type == "dna"){
                
                dt <- data.table::copy(read.genome(file = file, format = "fasta", ...))
        
                # here will come -> CDS prediction
                # then CDS -> protein translation
        }

        # WE COULD ALSO INSERT A IF-STATEMENT CHECKING THAT
        # IN CASE THE USER INSERTS AN AA-FASTA FILE,
        # THE TRANSLATION PROCESS IS BEING OMITTED


        # makedb
        dbname <- vector(mode = "character", length = 1)
        filename <- unlist(strsplit(file, .Platform$file.sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
        filename <- filename[length(filename)]


     if(makedb){
        if(!file.exists(file.path(tempdir(),"_blast_db"))){ 
                
                dir.create(file.path(tempdir(),"_blast_db")) 
                
        }
        
        currwd <- getwd()
        setwd(file.path(tempdir(),"_blast_db"))
        
        dbname <- paste0("blastdb_",filename,"_protein.fasta")
        
        seqinr::write.fasta( sequences = as.list(dt[ , aa]),
                             names     = dt[ , geneids],
                             nbchar    = 80, 
                             open      = "w",
                             file.out  = dbname )
        
        tryCatch({
                
                if(is.null(path)){
                        system(paste0("makeblastdb -in ", dbname,
                               " -input_type fasta -dbtype ",db_type," -hash_index"))
                        
                } else {
                        system(paste0("export PATH=",path,"; makeblastdb -in ",
                               dbname," -input_type fasta -dbtype ",db_type," -hash_index"))
                }
                
                }, error = function(e){ stop("makeblastdb did not work properly. The default parameters are: ","\n",
                                   "-input_type fasta -dbtype prot .","\n","Please check that you really want to work with a protein database.","\n",
                                   "Additionally check: ",dbname," .")}
        )
        
        # return to global working directory
        setwd(file.path(currwd))
     }

     return(list(stats::na.omit(dt), dbname))
}


