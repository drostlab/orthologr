#' @title Interface function to ProteinOrtho
#' @description This function takes a query organism and a set of subject organisms
#' and performs orthology inference using ProteinOrtho as orthology detection program (methodology).
#' @param query_file a character string specifying the path to the sequence file of interest (query organism).
#' @param subject_files a character string specifying the paths to the sequence files of interest (subject organisms).
#' @param po_params a character string specifying additional parameters that shall be handed to the PorteoOrtho call,
#' e.g. \code{po_params} = \code{"-identity=35 -cov=75"}. See \url{https://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html}
#' for parameter details.
#' @param eval a numeric value specifying the E-Value cutoff for BLAST hit detection.
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are: "protein", or "dna". ProteinOrtho can only handle protein or nucleotide sequences stored in fasta files.
#' @param format when using ProteinOrtho you can only specify \code{format} = \code{"fasta"}.
#' @param po_path a character string specifying the execution path to the \code{proteinortho5.pl} file (in case it is not sored
#' in the standard execution path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' used to run PrteinOrtho.
#' @param delete_files a boolean value specifying whether the folder '_ProteinOrtho' that stored the
#' ProteinOrtho output files shall be removed after the analysis. Default is \code{delete_files} = \code{FALSE}.
#' @details This function provides an interface to the ProteinOrtho program
#' and performs orthology inference using the ProteinOrtho methodology.
#' 
#' ProteinOrtho is a orthology inference program to detect orthologous genes within different species. 
#' For this purpose it compares similarities of given gene sequences and clusters them to find significant groups. 
#' To enhance the prediction accuracy, the relative order of genes (synteny) can be used as additional feature for the discrimination of orthologs.
#' (resource: \url{https://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html}).
#' 
#' @author Hajk-Georg Drost
#' @references
#' 
#' Lechner M et al. (2011) Proteinortho: Detection of (Co-)orthologs in large-scale analysis. BMC Bioinformatics 2011, 12:124.
#' 
#' \url{https://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html}
#' 
#' @examples \dontrun{
#' 
#' # finding orthologs between Arabidopsis thaliana and Arabidopsis lyrata genes
#' ProteinOrtho(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'              subject_files = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'              seq_type = "protein",format = "fasta")
#'              
#' # it is also possible to run ProteinOrtho using the synteny option
#' # make sure you can provide a *.gff file
#' # see https://www.bioinf.uni-leipzig.de/Software/proteinortho/manual.html for details
#' ProteinOrtho(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'              subject_files = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'              seq_type = "protein",format = "fasta",po_params = "-synteny")
#'    
#'              
#' # in case you need to specify the BLAST path or want to add 
#' # additional parameters to the BLAST run, please use the 'po_params' argument
#' ProteinOrtho(query_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
#'              subject_files = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
#'              seq_type = "protein",format = "fasta",
#'              po_params  "-blastpath='path/to/blastp' -blastParameters='seq no'")
#'              
#'              
#'                                  
#' }
#' @export

ProteinOrtho <- function(query_file, subject_files, po_params = NULL,eval = "1E-5",
                         seq_type = "protein", format = "fasta", po_path = NULL, 
                         comp_cores = 1, delete_files = FALSE){
        
        
        if(format != "fasta")
                stop("ProteinOrtho only supports fasta files storing protein sequences *.faa or nucleotide sequences *.fna")
        
        if(!is.element(seq_type,c("protein","dna")))
                stop("ProteinOrtho only supports protein sequences or nucleotide sequences. Please choose: 'protein' or 'dna'.")
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!file.exists(paste0("_ProteinOrtho",f_sep))){
                
                dir.create("_ProteinOrtho")
        }
        
        currwd <- getwd()
        setwd(file.path(currwd, "_ProteinOrtho"))
        
        # determine the number of cores on a multicore machine
        cores <- parallel::detectCores()
        
        # in case one tries to use more cores than are available
        if(comp_cores > cores)
                stop("You chose more cores than are available on your machine.")
        
        if(is.null(po_path)){
                
                
                if(is.null(po_params)){
                        
                        system(paste0("proteinortho5.pl -cpus=",comp_cores,
                                      " -project=ProteinOrtho ","-e=",eval," ",query_file," ",
                                      paste0(subject_files," ")))
                }
                
                if(!is.null(po_params)){
                        
                        system(paste0("proteinortho5.pl -cpus=",comp_cores,
                                      " -project=ProteinOrtho ","-e=",eval," ",po_params," ",
                                      query_file," ",paste0(subject_files," ")))
                }
        
        } else {
                
                
                if(is.null(po_params)){
                        
                        system(paste0(po_path,f_sep,"proteinortho5.pl -cpus=",comp_cores,
                                      " -project=ProteinOrtho ","-e=",eval," ",query_file," ",
                                      paste0(subject_files," ")))
                }
                
                if(!is.null(po_params)){
                        
                        system(paste0(po_path,f_sep,"proteinortho5.pl -cpus=",comp_cores,
                                      " -project=ProteinOrtho ","-e=",eval," ",po_params," ",
                                      query_file," ",paste0(subject_files," ")))
                }
                
                
        }
        
        setwd(currwd)
        
        tryCatch(
                {
        # read the header of the ProteinOrtho output file
        PO_tbl_header <- strsplit(readLines(paste0("_ProteinOrtho",f_sep,"ProteinOrtho.blast-graph"),n=3),"\t")
        PO_tbl_header <-  lapply(PO_tbl_header,stringr::str_replace_all,"#","")
        
        # store the output table of ProteinOrtho
        ProteinOrtho_tbl <- data.table::fread(paste0("_ProteinOrtho",f_sep,"ProteinOrtho.blast-graph"),sep = "\t",skip = 3)
        data.table::setnames(ProteinOrtho_tbl,old = paste0("V",1:dim(ProteinOrtho_tbl)[2]),
                             new = c(PO_tbl_header[[3]],PO_tbl_header[[2]][-c(1:length(PO_tbl_header[[3]]))]))
        data.table::setkeyv(ProteinOrtho_tbl,PO_tbl_header[[3]])
        
        if(delete_files)
                unlink("_ProteinOrtho",recursive = TRUE, force = TRUE)
        
        
        return(ProteinOrtho_tbl)
        
                }, error = function(e) stop(paste("The ProteinOrtho interface call did not terminate properly.",
                                            "Please make sure you passed all parameters correctly to ProteinOrtho.",
                                            "Did you provide a *.gff file in case you used the synteny option?",sep="\n"))
        )
        

}









