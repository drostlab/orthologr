#' @title Function to get a codon alignment.
#' @description This function takes a protein alignment and the corresponding coding sequences 
#' both given as files and uses the method specified by tool to create the corresponding codon
#' alignment.
#' @param file_aln a character string specifying the path to the file storing the protein alignment in CLUSTAL or FASTA format.
#' @param file_nuc a character string specifying the path to the file storing the coding sequences in multiple FASTA format.
#' @param format a character string specifying the file format used to store the codon alignment, e.g. "fasta", "clustal".
#' @param tool a character string specifying the program that should be used e.g. "pal2nal". 
#' @param codon_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{codon_aln_name} = \code{NULL} denoting a default name: 'tool_name.aln' .
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common codon alignment tools such as "PAL2NAL".
#' @examples \dontrun{
#' # performing a codon alignment using PAL2NAL
#' codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#'                        file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'), 
#'                        format = "clustal", tool = "pal2nal", get_aln = TRUE)
#' }
#' @return if get_aln is TRUE an object of class alignment of the seqinr package.
#' @export
codon_aln <- function(file_aln, file_nuc, format = "clustal", 
                      tool = "pal2nal", codon_aln_name = NULL, 
                      get_aln = FALSE, quiet = FALSE){
        
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        # the Pal2Nal program is stored as executable within
        # the R package environment: 'exec' folder
        # this is not an elegant version, only motivated by this discussion:
        # http://stackoverflow.com/questions/13463103/inst-and-extdata-folders-in-r-packaging
        # is there a better choice?
        path <- system.file(paste0("pal2nal",f_sep,"pal2nal.v14",f_sep,"pal2nal.pl"), package = "orthologr", mustWork = TRUE)
        #path <- "/exec/pal2nal.v14/"
        
        if(!is.element(tool,c("pal2nal")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(format,c("clustal", "fasta")))
                stop("Please choose a format that is supported by this function.")
        
        if(!file.exists(paste0("_alignment",f_sep))){
                
                dir.create("_alignment")
        }
        
        
        if(!file.exists(paste0("_alignment",f_sep,"codon_aln",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"codon_aln"))
        }
        
        if(is.null(codon_aln_name))
                file.out <- paste0("_alignment",f_sep,"codon_aln",f_sep,tool,".aln")
        
        if(!is.null(codon_aln_name))
                file.out <- paste0("_alignment",f_sep,"codon_aln",f_sep,codon_aln_name,"_",tool,".aln")
        
        # RIGHT NOW EACH NEW RUN OF THE FUNCTION OVERWRITES
        # THE EXISTING *.aln FILE
        # IN FUTURE VERSIONS WE SHOULD TRY TO KEEP (STORE) EXISTING FILES
        # AND RENAME THE NEW ONES -> session keys?
        
        # test whether the connection to pal2nal works
    
                if(tool == "pal2nal"){
                
                        tryCatch(
                        {   
                                system(paste0(path," ",file_aln," ",file_nuc," -output ",format," >",file.out))
                
                        },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                           "... the interface call did not work properly.") ) }
                        
                        )
                }

        if(!quiet){print(paste0("Codon Alignment successfully written in ",file.out,"."))}
        
                
        if(get_aln){
                
                tryCatch(
                        {
                                dna_aln <- seqinr::read.alignment(file = file.out, format = "clustal")
                                return(dna_aln)
               
                        }, error = function(){print(paste0("Something went wront with Pal2Nal.pl .\n",
                                                           file.out, " could not be read properly."))}
                )
        }

}

