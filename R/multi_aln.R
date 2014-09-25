#' @title An interface function to common multiple alignment tools.
#' @description This function takes a multiple FASTA file containing DNA or amino acid sequences
#' that shall be aligned and computes a multiple alignment using a defined multiple alignment tool.
#' @param file a character string specifying the path to the file storing the sequences in FASTA format.
#' @param tool a character string specifying the program that should be used: "clustalw", "t_coffee", "muscle", "clustalo", and "mafft". 
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param multi_aln_name a character string specifying the name of the stored alignment file. 
#' Default is \code{multi_aln_name} = \code{NULL} denoting a default name: 'tool_name.aln' .
#' @param clustalw.params a character string listing the input paramters that shall be passed to the executing clustalw program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running clustalw. Example: clustalw.params = "-PWMATRIX=BLOSUM -TYPE=PROTEIN".
#' @param t_coffee.params a character string listing the input paramters that shall be passed to the executing t_coffee program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running t_coffee. Here the default '-mode' is 'accurate' (only works with protein sequences). Example: t_coffee.params = "-mode expresso".
#' @param muscle.params a character string listing the input paramters that shall be passed to the executing muscle program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running muscle. Example: muscle.params = "-diags -clwstrict".
#' @param clustalo.params a character string listing the input paramters that shall be passed to the executing clustalo program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running clustalo. Example: clustalo.params = "--threads 2 --outfmt clu"
#' @param mafft.params a character string listing the input paramters that shall be passed to the executing mafft program. Default is \code{NULL}, implicating
#' that a set of default parameters is used when running mafft. Example: mafft.params = "--maxiterate 1 --thread 2 --clustalout".
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function provides an interface between R and common multiple alignment programs
#' such as "clustalw", "t_coffee", "muscle", "clustalo", and "mafft".
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
#' T_COFFEE :
#' 
#' In case you use the default path to the t_coffee program,
#' the following calls to clustalw should work properly on your system:
#' 
#' system("t_coffee -version")
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 't_coffee' execution path on your system:
#' 
#' system("path/to/t_coffee/t_coffee -version")
#' 
#' 
#' MUSCLE : 
#' 
#' In case you use the default path to the muscle program,
#' the following calls to muscle should work properly on your system:
#' 
#' system("muscle -help")
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'muscle' execution path on your system:
#' 
#' system("path/to/muscle/muscle -help")
#' 
#' 
#' #' CLUSTALO : 
#' 
#' In case you use the default path to the clustalo program,
#' the following calls to clustalo should work properly on your system:
#' 
#' system("clustalo --help")
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'clustalo' execution path on your system:
#' 
#' system("path/to/clustalo/clustalo --help")
#' 
#' 
#' MAFFT : 
#' 
#' In case you use the default path to the mafft program,
#' the following calls to mafft should work properly on your system:
#' 
#' system("mafft -help")
#' 
#' In case this procedures doesn't work properly, please use the \code{path} argument
#' to specify the 'mafft' execution path on your system:
#' 
#' system("path/to/mafft/mafft -help")
#' 
#' @note Note that when using the 'clustalw.params', ... parameters, make sure the corresponding alignment tool
#' returns a file in clustal format *.aln. This is only important when \code{get_aln} = \code{TRUE}.
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
#' # running clustalw using additional parameters
#' # details: system("clustalw2 -help")
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "clustalw", get_aln = TRUE, 
#'           clustalw.params = "-PWMATRIX=BLOSUM -TYPE=PROTEIN")
#'           
#'           
#' # T_COFFEE Example:
#' 
#' # in case the default execution path of t_coffee runs properly on your system
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "t_coffee", get_aln = TRUE)
#' 
#' # in case the default execution path of t_coffee is not set within the default path
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "t_coffee", get_aln = TRUE, path = "path/to/t_coffee/")
#' 
#' # running t_coffee using additional parameters
#' # details: http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "t_coffee", get_aln = TRUE,
#'           t_coffee.params = "-mode expresso")
#'           
#'
#'
#' # MUSCLE Example:  
#' 
#' # in case the default execution path of muscle runs properly on your system
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "muscle", get_aln = TRUE)
#' 
#' # in case the default execution path of muscle is not set within the default path
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "muscle", get_aln = TRUE, path = "path/to/muscle/")
#' 
#' # running muscle using additional parameters
#' # details: http://www.drive5.com/muscle/manual/
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "muscle", get_aln = TRUE,
#'           muscle.params = "-diags -clwstrict") 
#'           
#' 
#' # CLUSTALO Example:  
#' 
#' # in case the default execution path of clustalo runs properly on your system
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "clustalo", get_aln = TRUE)
#' 
#' # in case the default execution path of clustalo is not set within the default path
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "clustalo", get_aln = TRUE, path = "path/to/clustalo/")
#' 
#' # running clustalo using additional parameters
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "clustalo", get_aln = TRUE,
#'           clustalo.params = "--outfmt clu")         
#'           
#'           
#'                                             
#' # MAFFT Example:  
#' 
#' # in case the default execution path of mafft runs properly on your system
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "mafft", get_aln = TRUE)
#' 
#' # in case the default execution path of mafft is not set within the default path
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "mafft", get_aln = TRUE, path = "path/to/mafft/")
#' 
#' # running mafft using additional parameters
#' # details: http://www.drive5.com/mafft/manual/
#' multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#'           tool = "mafft", get_aln = TRUE,
#'           mafft.params = "--maxiterate 1 --clustalout")         
#'           
#'                                                                                   
#' }
#' @references 
#' 
#' CLUSTALW:
#' 
#' Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F,
#' Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ, Higgins DG. (2007).
#' Clustal W and Clustal X version 2.0. Bioinformatics, 23, 2947-2948.
#' 
#' http://www.clustal.org/clustal2/
#' 
#' http://www.ebi.ac.uk/Tools/msa/clustalw2/help/
#' 
#' 
#' T_COFFEE
#' 
#' T-Coffee: A novel method for multiple sequence alignments. 
#' Notredame, Higgins, Heringa, JMB, 302(205-217). 2000.
#' 
#' http://www.tcoffee.org/Projects/tcoffee/
#' 
#' http://www.tcoffee.org/Projects/tcoffee/documentation/t_coffee_tutorial.pdf
#' 
#' 
#' MUSCLE:
#' 
#' Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res. 32(5):1792-1797. 
#' 
#' Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced time and space complexity. BMC Bioinformatics, (5) 113. 
#' 
#' http://www.drive5.com/muscle/
#' 
#' http://www.drive5.com/muscle/manual/
#' 
#' 
#' CLUSTALO:
#' 
#' Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). 
#' Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75
#' 
#' http://www.clustal.org/omega/
#' 
#' http://www.clustal.org/omega/README
#' 
#' MAFFT : 
#' 
#' Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)  
#' MAFFT multiple sequence alignment software version 7: improvements in performance and usability. 
#' 
#' http://mafft.cbrc.jp/alignment/software/
#' 
#' http://mafft.cbrc.jp/alignment/software/manual/manual.html
#' 
#' http://mafft.cbrc.jp/alignment/software/tips0.html
#' 
#' @return In case the argument \code{get_aln} is set \code{TRUE}, an object of class alignment of the seqinr package is returned.
#' @export
multi_aln <- function(file, tool, get_aln = FALSE, path = NULL,
                      multi_aln_name = NULL, clustalw.params = NULL, 
                      t_coffee.params = NULL, muscle.params = NULL, 
                      clustalo.params = NULL, mafft.params = NULL, 
                      quiet = FALSE){
        
        if(!is.element(tool,c("clustalw", "t_coffee", "muscle", "clustalo","mafft")))
                stop("Please choose a tool that is supported by this function.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!file.exists(paste0("_alignment",f_sep))){
                
                dir.create("_alignment")
        }
        
        if(!file.exists(paste0("_alignment",f_sep,"multi_aln",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"multi_aln"))
        }
        
        if(is.null(multi_aln_name)){
                
                file.out <- paste0("_alignment",f_sep,"multi_aln",f_sep,tool,".aln")
        }

        if(!is.null(multi_aln_name)){
         
                file.out <- paste0("_alignment",f_sep,"multi_aln",f_sep,multi_aln_name,"_",tool,".aln")    
        }
        
# does not work as expected, could be included in another way or not.
#         if(quiet){
#                 clustalw.params = paste0(clustalw.params, " -QUIET")
#                 mafft.params = paste0(mafft.params, " --quiet")
#                 muscle.params = paste0(muscle.params, " -quiet")
#                 # clustalo - no param found to limit the output
#         }

        
        # RIGHT NOW EACH NEW RUN OF THE FUNCTION OVERWRITES
        # THE EXISTING *.aln FILE
        # IN FUTURE VERSIONS WE SHOULD TRY TO KEEP (STORE) EXISTING FILES
        # AND RENAME THE NEW ONES
        
        if(tool == "clustalw"){
                
                # find out on what kind of OS ClustalW is running
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
                               
                               # if no specific parameters are set,
                               # then use the default parameters
                               if(is.null(clustalw.params)){
                                       
                                        # right now only the default parameters for args: "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN", "-GAPEXT" are being used
                                        system(paste0(call_clustalw," -infile=",file," -outfile=",file.out," -quiet"))
                               } else {
                                       # add additional parameters when running clustalw
                                        system(paste0(call_clustalw," -infile=",file," -outfile=",file.out," ",clustalw.params))                                       
                               }
                        } else {
                                
                                # if no specific parameters are set,
                                # then use the default parameters
                                if(is.null(clustalw.params)){
                                        
                                        # use the default parameters when running clustalw
                                        system(
                                                paste0("export PATH=$PATH:",path,"; ",call_clustalw," -infile=",
                                                file," -outfile=",file.out," -quiet")
                                        )
                                } else {
                                        
                                        # add additional parameters when running clustalw
                                        system(
                                                paste0("export PATH=$PATH:",path,"; ",call_clustalw," -infile=",
                                                       file," -outfile=",file.out," ",clustalw.params)
                                        )    
                                }
                        }
                },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                   "... the interface call did not work properly.") )}
                )
                
                if(! quiet){print(paste0("Multiple Alignment successfully written in ",file.out,"."))}
                
                
                if(get_aln){
                        
                        tryCatch(
                                {
                                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                                        return(aln)
                                        
                                }, error = function(){print(paste0("Something went wront with ",tool," .\n",
                                                                   file.out, " could not be read properly."))}
                        )
                }
        }
 
        if(tool == "t_coffee"){
                
                # test whether the connection to t_coffee works
                tryCatch(
                {
                        if(is.null(path)){
                                
                                if(is.null(t_coffee.params)){
                                        
                                        # use the default parameters when running t_coffee
                                        # perform an accurate alignment which is very accurate, but slow
                                        # http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
                                        system(paste0("t_coffee -infile ",file," -mode accurate"," -outfile ",file.out))
                                } else {
                                        
                                        # add additional parameters when running t_coffee
                                        system(paste0("t_coffee -infile ",file," ",t_coffee.params," -outfile ",file.out))
                                        
                                }
                        } else {
                                
                                if(is.null(t_coffee.params)){
                                        
                                        # use the default parameters when running t_coffee
                                        # perform a accurate alignment which is very accurate, but slow
                                        # http://www.tcoffee.org/Projects/tcoffee/#DOCUMENTATION
                                        system(paste0("export PATH=$PATH:",path,"; ","t_coffee -infile ",
                                                 file," -mode accurate"," -outfile ",file.out))
                                } else {
                                        
                                        # add additional parameters when running t_coffee
                                        system(paste0("export PATH=$PATH:",path,"; ","t_coffee -infile ",
                                                      file," ",t_coffee.params," -outfile ",file.out))
                                }
                        }
                
                },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                   "... the interface call did not work properly.") )}
                )
                
                if(! quiet){print(paste0("Multiple Alignment successfully written in ",file.out,"."))}
                
                
                if(get_aln){
                        
                        tryCatch(
                                {
                                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                                        return(aln)
                                }, error = function(){print(paste0("Something went wront with ",tool," .\n",
                                                                   file.out, " could not be read properly."))}
                        )
                }
        }
 
         if(tool == "muscle"){
                 
                 # test whether the connection to muscle works
                 tryCatch(
                 {        
                          if(is.null(path)){
                                  if(is.null(muscle.params)){
                                          
                                          # use the default parameters when running muscle
                                          # write output into clustalw format using the -clwstrict argument
                                           system(paste0("muscle -in ",file," -out ",file.out, " -clwstrict"," -quiet"))
                                           
                                  } else {
                                          
                                          # add additional parameters when running muscle
                                          # write output into clustalw format using the -clwstrict argument
                                           system(paste0("muscle -in ",file," ",muscle.params," -out ",file.out))
                                           
                                  }
                 
                           } else {
                                   if(is.null(muscle.params)){
                                           
                                           # use the default parameters when running muscle
                                           # write output into clustalw format using the -clwstrict argument
                                           system(paste0("export PATH=$PATH:",path,"; ","muscle -in ",
                                                  file," -out ",file.out," -clwstrict"," -quiet"))
                                           
                                   } else {
                                           
                                           # add additional parameters when running muscle
                                           # write output into clustalw format using the -clwstrict argument
                                           system(paste0("export PATH=$PATH:",path,"; ","muscle -in ",
                                                         file," ",muscle.params," -out ",file.out))
                                           
                                   }
                         
                           }
                 
                 },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                    "... the interface call did not work properly.") )}
                 )
                 
                 if(! quiet){print(paste0("Multiple Alignment successfully written in ",file.out,"."))}
                 
                 
                 if(get_aln){
                         
                         tryCatch(
                                 {
                                         aln <- seqinr::read.alignment(file.out, format = "clustal")
                                         return(aln)
                                         
                                 }, error = function(){print(paste0("Something went wront with ",tool," .\n",
                                                                     file.out, " could not be read properly."))}
                         )
                 }
         }
 
         if(tool == "clustalo"){
                 
                 # test whether the connection to clustalo works
                 tryCatch(
                 {       
                         if(is.null(path)){
                                 
                                 if(is.null(clustalo.params)){
                                         
                                         # use the default parameters when running clustalo
                                         system(paste0("clustalo -i ",file," -o ",file.out," --outfmt clu --force"))
                                 
                                 } else {
                                         
                                         # add additional parameters when running clustalo
                                         system(paste0("clustalo -i ",file," -o ",file.out," --force ",clustalo.params))
                                 }
                 
                         } else {
                                 if(is.null(clustalo.params)){
                                         
                                         # use the default parameters when running clustalo
                                         system(paste0("export PATH=$PATH:",path,"; ","clustalo -i ",
                                                file," -o ",file.out," --outfmt clustal --force"))
                                 
                                 } else {
                                         
                                         # add additional parameters when running clustalo
                                         system(paste0("export PATH=$PATH:",path,"; ","clustalo -i ",
                                                       file," -o ",file.out," --force ",clustalo.params))
                                         
                                 }
                         
                        }
                 
                 },error = function(){ print(paste0("Please check the correct path to ",tool,
                                                    "... the interface call did not work properly.") )}
                 )
                 
                 if(! quiet){print(paste0("Multiple Alignment successfully written in ",file.out,"."))}
                 
                         
                 
                 if(get_aln){
                         
                         tryCatch(
                                 {
                                         aln <- seqinr::read.alignment(file.out, format = "clustal")
                                         return(aln)
                                         
                         }, error = function(){print(paste0("Something went wront with ",tool," .\n",
                                                                       file.out, " could not be read properly."))}
                         )
                 }
         }
        
        if(tool == "mafft"){
                
                # test whether the connection to mafft works
                tryCatch(
                {
                 if(is.null(path)){
                
                         # if no specific parameters are set,
                         # then use the default parameters
                         if(is.null(mafft.params)){
                                
                                # use the default parameters when running mafft
                                system(paste0("mafft --quiet --clustalout ",file," >",file.out))
                                
                         } else {
                                
                                # add additional parameters when running mafft
                                system(paste0("mafft ",mafft.params," ",file," >",file.out))                                       
                         }
                } else {
                
                         # if no specific parameters are set,
                         # then use the default parameters
                         if(is.null(mafft.params)){
                                
                                # use the default parameters when running mafft
                                system(
                                        paste0("export PATH=$PATH:",path,"; ","mafft --quiet --clustalout ",
                                               file," >",file.out)
                                       )
                         } else {
                                
                                # add additional parameters when running mafft
                                system(
                                        paste0("export PATH=$PATH:",path,"; ","mafft"," ",
                                               mafft.params," ",file," >",file.out)
                                )    
                        }
                }
               },error = function(){ print(paste0("Please check the correct path to ",tool,
                                   "... the interface call did not work properly.") )}
                )
               
               if(! quiet){print(paste0("Multiple Alignment successfully written in ",file.out,"."))}
               
               
               ## LINUX ERROR 
               #                  There is a problem in the configuration of your shell.
               #                  Check the MAFFT_BINARIES environmental variable by
               #                  $ echo $MAFFT_BINARIES
               #                  This variable must be *unset*, unless you have installed MAFFT
               #                  with a special configuration.  To unset this variable, type
               #                  $ unset MAFFT_BINARIES or $ unsetenv MAFFT_BINARIES                  
               #                  To keep this change permanently, edit setting files
               #                  (.bash_profile, .profile, .cshrc, etc) in your home directory
               #                  to delete the MAFFT_BINARIES line.
               

                if(get_aln){
                        
                        tryCatch(
                                {
                                        aln <- seqinr::read.alignment(file.out, format = "clustal")
                                        return(aln)
                                }, error = function(){print(paste0("Something went wront with ",tool," .\n",
                                                                   file.out, " could not be read properly."))}
                        )
                }
        }

}















