#' @title An interface function to common multiple alignment tools.
#' @description This function takes a multiple FASTA file containing DNA or amino acid sequences
#' that shall be aligned and computes a multiple alignment using a defined multiple alignment tool.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program that should be used: "clustalw", "tcoffee", "muscle", and "clustalo". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between R and common multiple alignment programs
#' such as "clustalw", "tcoffee", "muscle", and "clustalo".
#' 
#' CLUSTALW : 
#' 
#' Different operating systems perform different execution calls to the clustalw program:
#' 
#' MacOS: 'clustalw2', Linux: 'clustalw', Windows: 'clustalw2.exe'
#' 
#' In case you use the default path to the clustalw program, depending on your operating system,
#' the following calls to clustalw should work properly on your system:
#' 
#' MacOS: system("clustalw2 -help"), Linux: system("clustalw -help"), Windows: system("clustalw2.exe -help")
#' 
#' In case these procedures don't work properly, please use the \code{path} argument
#' to specify the 'clustalw' execution path on your system:
#' 
#' 
#' MacOS: system("path/to/clustalw/clustalw2 -help"), Linux: system("path/to/clustalw/clustalw -help"), Windows: system("path/to/clustalw/clustalw2.exe -help")
#' 
#' 
#' @examples \dontrun{
#' 
#' # CLUSTALW Example:
#' 
#' # in case the default execution path of clustalw runs properly on your system
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "clustalw", get_aln = TRUE)
#' 
#' # in case the default execution path of clustalw is not set within the default path
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'), 
#'           tool = "clustalw", get_aln = TRUE, path = "path/to/clustalw/")
#' 
#' 
#' }
#' @return If get_aln is TRUE an object of class alignment of the seqinr package.
#' @export
multi_aln <- function(file, tool, get_aln = "FALSE", path = NULL){
        
        if(!is.element(tool,c("clustalw", "tcoffee", "muscle", "clustalo")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!file.exists("_alignment/")){
                
                dir.create("_alignment")
        }
        
        file.out <- paste0("_alignment/",tool,".aln")
        
        # RIGHT NOW EACH NEW RUN OF THE FUNCTION OVERWRITES
        # THE EXISTING *.aln FILE
        # IN FUTURE VERSIONS WE SHOULD TRY TO KEEP (STORE) EXISTING FILES
        # AND RENAME THE NEW ONES
        
        if(tool == "clustalw"){
                
                
                # find out on what kind of OS this running
                operating_sys <- Sys.info()[1]
                
                if (operating_sys == "Darwin") 
                        call_clustalw <- "clustalw2"
                
                if (operating_sys == "Linux")
                        call_clustalw <- "clustalw"
                
                if (operating_sys == "Windows") 
                        call_clustalw <- "clustalw2.exe"
                
                
                # test whether the connection to clustalw works
                tryCatch(
                {
                       if(is.null(path)){
                               
                               # right now only the default parameters for args: "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN", "-GAPEXT" are used
                                system(paste0(call_clustalw," -infile=",file," -outfile=",file.out," -quiet"))
                        
                        } else {
                                # right now only the default parameters for args: "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN", "-GAPEXT" are used
                                system(
                                        paste0("export PATH=$PATH:",path,"; ",call_clustalw," -infile=",
                                               file," -outfile=",file.out," -quiet")
                                )
                        }
                },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                   "... the interface call did not work properly.") )}
                , finally = print(paste0("Multiple Alignment successfully written in ",file.out,"."))
                
                )
                
                if(get_aln){
                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                        return(aln)
                }
        }
 
        if(tool == "tcoffee"){
                
                # test whether the connection to t_coffee works
                tryCatch(
                {
                        if(is.null(path)){
                        
                                system(paste0("t_coffee ",file," >",file.out))
                        
                        } else {
                        
                                system(paste0("export PATH=$PATH:",path,"; ","t_coffee ",
                                              file," >",file.out))
                        
                        }
                
                },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                   "... the interface call did not work properly.") )}
                , finally = print(paste0("Multiple Alignment successfully written in ",file.out,"."))
                
                )
                
                if(get_aln){
                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                        return(aln)
                }
        }
 
         if(tool == "muscle"){
                 
                 # test whether the connection to muscle works
                 tryCatch(
                 {        
                          if(is.null(path)){
                         
                                   system(paste0("muscle -in ",file," -out ",file.out))
                 
                           } else {
                         
                                   system(paste0("export PATH=$PATH:",path,"; ","muscle -in ",
                                       file," -out ",file.out))
                         
                           }
                 
                 },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                    "... the interface call did not work properly.") )}
                 , finally = print(paste0("Multiple Alignment successfully written in ",file.out,"."))
                 
                 )
                 
                 if(get_aln){
                         aln <- seqinr::read.alignment(file.out, format = "fasta")
                         return(aln)
                 }
         }
 
         if(tool == "clustalo"){
                 
                 # test whether the connection to clustalo works
                 tryCatch(
                 {       
                         if(is.null(path)){
                         
                                 system(paste0("clustalo -i ",file," -o ",file.out," --outfmt clustal --force"))
                 
                         } else {
                         
                                 system(paste0("export PATH=$PATH:",path,"; ","clustalo -i ",
                                               file," -o ",file.out," --outfmt clustal"))
                         
                        }
                 
                 },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                    "... the interface call did not work properly.") )}
                 , finally = print(paste0("Multiple Alignment successfully written in ",file.out,"."))
                 
                 )
         
                 
                 if(get_aln){
                         aln <- seqinr::read.alignment(file.out, format = "clustal")
                         return(aln)
                 }
         }

}




