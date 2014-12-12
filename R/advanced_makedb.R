#' @title Create a BLASTable database with makeblastdb (advanced options)
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



