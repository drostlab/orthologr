#' @title Function to calculate the synonymous vs nonsynonymous substitutionrates for two organisms.
#' @description This function takes the CDS files of two organisms of interest (query_file and subject_file)
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms. 
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param seq_type a character string specifying the sequence type stored in the input file.
#' Options are are: "cds", "protein", or "dna". In case of "cds", sequence are translated to protein sequences,
#' in case of "dna", cds prediction is performed on the corresponding sequences which subsequently are
#' translated to protein sequences. Default is \code{seq_type} = "cds".
#' @param format a character string specifying the file format of the sequence file, e.g. "fasta", "gbk". See \code{\link{read.cds}},
#' \code{\link{read.genome}}, \code{\link{read.proteome}} for more details.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed
#' to detect orthologous genes. Default is \code{ortho_detection} = "RBH" (BLAST reciprocal best hit).
#' Further methods are: "BH" (BLAST best hit), "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), "OrthoMCL, "IP" (InParanoid).
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type} = "multiple" or \code{aa_aln_type} = "pairwise".
#' Default is \code{aa_aln_type} = "multiple".
#' @param aa_aln_tool a character string specifying the program that should be used e.g. "clustalw".
#' @param aa_aln_path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param aa_aln_params  a character string specifying additional parameters that shall be passed to the selected alignment tool. Default is \code{aa_aln_params} = \code{NULL} 
#' (no addintional parameters are passed to the selected alignment tool).
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used. Default is \code{codon_aln_tool} = \code{"pal2nal"}.
#' Right now only "pal2nal" can be selected as codon alignment tool.
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine. 
#' @param quiet a logical value specifying whether the output of the corresponding alignment tool shall be printed out to the console.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details 
#' 
#' #' The dNdS estimation methods available in this function are:
#' 
#' - "Li" : Li's method (1993) -> provided by the ape package
#' 
#' - "Comeron" : Comeron's method (1995)
#' 
#' dNdS estimation methods provided by KaKs_Calculator 1.2 :
#' 
#' Approximate Methods:
#' 
#' "NG": Nei, M. and Gojobori, T. (1986)
#' 
#' "LWL": Li, W.H., et al. (1985)
#' 
#' "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' 
#' "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' 
#' "YN": Yang, Z. and Nielsen, R. (2000)
#' 
#' "MYN" (Modified YN): Zhang, Z., et al. (2006)
#' 
#' Maximum-Likelihood Methods:
#' 
#' GY: Goldman, N. and Yang, Z. (1994)
#' 
#' MS (Model Selection), MA (Model Averaging): based on a set of candidate models defined by Posada, D. (2003) as follows.
#' 
#' MS (Model Selection) and MA (Model Averaging) over:
#' 
#' JC,
#' F81,
#' K2P, 
#' HKY,
#' TrNEF,
#' TrN,
#' K3P,
#' K3PUF,
#' TIMEF,
#' TIM,
#' TVMEF,
#' TVM,
#' SYM,
#' GTR
#' 
#' 
#' @return A data.table storing the dNdS values of the correspnding genes.
#' @examples \dontrun{
#' 
#' # get a dNdS table using:
#' # 1) reciprocal best hit for orthology inference (RBH)
#' # 2) clustalw for pairwise amino acid alignments
#' # 3) pal2nal for codon alignments
#' # 4) Yang, Z. and Nielsen, R. (2000) (YN) for dNdS estimation
#' # 5) single core processing 'comp_cores = 1'
#' dNdS(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#' subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#' ortho_detection = "RBH", aa_aln_type = "multiple",
#' aa_aln_tool = "clustalw", aa_aln_path = "/path/to/clustalw/",
#' codon_aln_tool = "pal2nal", dnds_est.method = "YN", comp_cores = 1)
#' 
#' # The same result can be obtained using multicore processing using: comp_cores = 2 or 3 or more ...
#' dNdS(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#' subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#' ortho_detection = "RBH", aa_aln_type = "multiple",
#' aa_aln_tool = "clustalw", aa_aln_path = "/path/to/clustalw/",
#' codon_aln_tool = "pal2nal", dnds_est.method = "YN", comp_cores = 2)
#' 
#' }
#' @seealso \code{\link{substitutionrate}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, \code{\link{blast_best}},
#' \code{\link{blast_rec}}, \code{\link{read.cds}}
#' @export
dNdS <- function(query_file, subject_file, seq_type = "protein",
                 format = "fasta", ortho_detection = "RBH", 
                 blast_path = NULL, aa_aln_type = "multiple", 
                 aa_aln_tool = "clustalw", aa_aln_path = NULL, 
                 aa_aln_params = NULL, codon_aln_tool = "pal2nal", 
                 dnds_est.method = "YN", comp_cores = 1, 
                 quiet = FALSE, clean_folders = FALSE){
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!is.ortho_detection_method(ortho_detection))
                stop("Please choose a orthology detection method that is supported by this function.")
        
        if(!is.element(aa_aln_type,c("multiple","pairwise")))
                stop("Please choose a supported alignement type: 'multiple' or 'pairwise'")
        
        if(aa_aln_type == "multiple"){
                if(!is.multiple_aln_tool(aa_aln_tool))
                        stop("Please choose a multiple alignment tool that is supported by this function or try to choose aa_aln_type = 'pairwise'.")
        }
        if(aa_aln_type == "pairwise"){
                if(!is.pairwise_aln_tool(aa_aln_tool))
                        stop("Please choose a pairwise alignment tool that is supported by this function or try to choose aa_aln_type = 'multiple'.")
        }

        if(!is.codon_aln_tool(codon_aln_tool))
                stop("Please choose a codon alignment tool that is supported by this function.")
        
        if(!is.dnds_est_method(dnds_est.method))
                stop("Please choose a dNdS estimation method that is supported by this function.")
        
        # blast each translated aminoacid sequence against the related database to get a 
        # hit table with pairs of geneids  
        
        if(ortho_detection == "BH"){
                
                # seq_type = "cds" -> dNdS() needs CDS files as input!
                hit.table <- data.table::copy( blast_best(query_file = query_file, subject_file = subject_file, 
                                   path = blast_path, comp_cores = comp_cores,
                                   seq_type = "cds", format = format))
                                                
                q_cds <- read.cds(file = query_file, format = format)
                s_cds <- read.cds(file = subject_file, format = format)
                
                q_aa <- read.proteome(file = paste0("_blast_db",f_sep,"blastinput.fasta"), format = "fasta")
                
                filename_subj <- unlist(strsplit(subject_file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename_subj <- filename_subj[length(filename_subj)]
                s_aa <- read.proteome(file = paste0("_blast_db",f_sep,"out_",filename_subj,"_translate.fasta"), format = "fasta")
    
        }
        
        if(ortho_detection == "RBH"){
                           
                # seq_type = "cds" -> dNdS() needs CDS files as input!
                hit.table <- data.table::copy(
                        blast_rec(query_file = query_file, subject_file = subject_file, 
                                   path = blast_path, comp_cores = comp_cores,
                                   seq_type = "cds", format = format))
                
                
                q_cds <- read.cds(file = query_file, format = format)
                s_cds <- read.cds(file = subject_file, format = format)
                
                
                filename_qry <- unlist(strsplit(query_file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename_qry <- filename_qry[length(filename_qry)]
                
                q_aa <- read.proteome(file = paste0("_blast_db",f_sep,"out_",filename_qry,"_translate.fasta"), format = "fasta")
                
                filename_subj <- unlist(strsplit(subject_file, f_sep, fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename_subj <- filename_subj[length(filename_subj)]
                
                s_aa <- read.proteome(file = paste0("_blast_db",f_sep,"out_",filename_subj,"_translate.fasta"), format = "fasta")
                
                
        }
    
        data.table::setnames(q_cds, old=c("geneids", "seqs"), new = c("query_id","query_cds"))
        data.table::setnames(s_cds, old=c("geneids", "seqs"), new = c("subject_id","subject_cds"))
        data.table::setnames(q_aa, old=c("geneids", "seqs"), new = c("query_id","query_aa"))
        data.table::setnames(s_aa, old=c("geneids", "seqs"), new = c("subject_id","subject_aa"))
        
        # joining all tables to a final table containing: query_id, subject_id, query_aa_seq, subject_aa_seq, query_cds_seq, and subject_cds_seq
        query_tbl <- dplyr::inner_join(dplyr::tbl_dt(q_cds), dplyr::tbl_dt(q_aa), by = "query_id")
        subject_tbl <- dplyr::inner_join(dplyr::tbl_dt(s_cds), dplyr::tbl_dt(s_aa), by = "subject_id")
        joint_query_tbl <- dplyr::inner_join(dplyr::tbl_dt(hit.table), dplyr::tbl_dt(query_tbl), by = "query_id")
        joint_subject_tbl <- dplyr::inner_join(dplyr::tbl_dt(hit.table), dplyr::tbl_dt(subject_tbl), by = "subject_id")
        final_tbl <- dplyr::inner_join(dplyr::tbl_dt(joint_query_tbl), dplyr::tbl_dt(joint_subject_tbl), by = c("query_id","subject_id"))
        
        if(nrow(final_tbl) == 0)
                stop("No orthologs could be found! Please check your input files!")
        
        
       # return the dNdS table for all query_ids and subject_ids
           return(compute_dnds(complete_tbl = final_tbl,
                       aa_aln_type = aa_aln_type,
                       aa_aln_tool = aa_aln_tool,
                       aa_aln_path = aa_aln_path,
                       codon_aln_tool = codon_aln_tool, 
                       dnds_est.method = dnds_est.method, quiet = quiet,
                       comp_cores = comp_cores, clean_folders = clean_folders)
                       )   
}




#' @title Function to calculate the synonymous vs nonsynonymous substitutionrate for a codon alignment.
#' @description This function takes a pairwise alignment as input file and estimates the
#' dNdS ratio of the corresponding alignment. Nevertheless, this function is a helper function for
#' \code{\link{dNdS}}. For dNdS computations you should use the function: \code{\link{dNdS}}.
#' @param file a character string specifying the path to a codon alignment file
#' @param est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li" .
#' Note, that when using "Comeron" as dNdS estimation method, the program 'gestimator' is used to compute the
#' corresponding dNdS values from a given alignment. The program 'gestimator' can only read "fasta" files,
#' hence it is important to use format = "fasta" when choosing est.method = "Comeron".
#' @param format a character string specifying the file format in which the alignment is stored:  
#' "mase", "clustal", "phylip", "fasta" , "msf"
#' @param quiet a logical value specifying whether the output of the coresponding interface shall be printed out.
#' @param kaks_calc.params a character string storing additional parameters for KaKs_Claculator 1.2 . Default is \code{NULL}. Example:
#' \code{kaks_calc.params} = "-m NG -m YN". 
#' @param subst_name a character string specifying the substitution name that shall be added to the internal folder path naming. 
#' Default is \code{subst_name} = \code{NULL}.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes a pairwise alignments file as input and estimates dNdS ratios
#' of the corresponding alignments using predefined dNdS estimation methods.
#' 
#' The dNdS estimation methods available in this function are:
#' 
#' - "Li" : Li's method (1993) -> provided by the ape package
#' 
#' - "Comeron" : Comeron's method (1995)
#' 
#' dNdS estimation methods provided by KaKs_Calculator 1.2 :
#' 
#' Approximate Methods:
#' 
#' "NG": Nei, M. and Gojobori, T. (1986)
#' 
#' "LWL": Li, W.H., et al. (1985)
#' 
#' "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
#' 
#' "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
#' 
#' "YN": Yang, Z. and Nielsen, R. (2000)
#' 
#' "MYN" (Modified YN): Zhang, Z., et al. (2006)
#' 
#' Maximum-Likelihood Methods:
#' 
#' GY: Goldman, N. and Yang, Z. (1994)
#' 
#' MS (Model Selection), MA (Model Averaging): based on a set of candidate models defined by Posada, D. (2003) as follows.
#' 
#' MS (Model Selection) and MA (Model Averaging) over:
#' 
#' JC,
#' F81,
#' K2P, 
#' HKY,
#' TrNEF,
#' TrN,
#' K3P,
#' K3PUF,
#' TIMEF,
#' TIM,
#' TVMEF,
#' TVM,
#' SYM,
#' GTR
#' 
#' @examples \dontrun{
#' 
#' # estimate the dNdS rate using Li's method
#' substitutionrate(system.file("seqs/aa_seqs.aln", package = "orthologr"),
#'                  est.method = "Li", format = "clustal")
#'  
#'  # estimate the dNdS rate using model averaging provided by the KaKs_Calculator 1.2 program
#'  substitutionrate(system.file("seqs/pal2nal.aln", package = "orthologr"), 
#'                   est.method = "MA", format = "fasta") 
#'                   
#'  # estimate the dNdS rate using Nei and Gojobori's method provided by the 
#'  # KaKs_Calculator 1.2 program
#'  substitutionrate(system.file("seqs/pal2nal.aln", package = "orthologr"), 
#'                   est.method = "NG", format = "fasta")     
#'   
#'   # estimate the dNdS rate using Nei and Gojobori's method AND Yang and Nielsen's 
#'   # method provided by the KaKs_Calculator 1.2 program 
#'   # for this purpose we choose: est.method = "kaks_calc" and kaks_calc.params = "-m NG -m YN"                 
#'  substitutionrate(system.file("seqs/pal2nal.aln", package = "orthologr"),
#'                   est.method = "kaks_calc", format = "fasta",kaks_calc.params = "-m NG -m YN")
#'                   
#'                                           
#' }
#' @references 
#' Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol., 36:96-99.
#' 
#' Charif, D. and Lobry, J.R. (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.
#' 
#' Thornton, K. (2003) libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19(17): 2325-2327
#' 
#' Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator:
#' Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263.
#' 
#' https://code.google.com/p/kaks-calculator/wiki/KaKs_Calculator
#' 
#' https://code.google.com/p/kaks-calculator/wiki/AXT
#' 
#' @return A data.table storing the query_id, subject_id, dN, dS, and dNdS values or 
#' a data.table storing the query_id, method, dN, dS, and dNdS values when using KaKs_Calculator.
#' @seealso \code{\link{dNdS}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, \code{\link{blast_best}},
#' \code{\link{blast_rec}}, \code{\link{read.cds}}
#' @export
substitutionrate <- function(file, est.method, format = "fasta",
                             quiet = FALSE, kaks_calc.params = NULL,
                             subst_name = NULL){
        
        # dNdS estimation methods provided by the KaKs_Calculator 1.2 program
        kaks_calc_methods <- c("MA","MS","NG","LWL","LPB","MLWL","YN","MYN","GY","kaks_calc")
        
        if(!is.dnds_est_method(est.method))
                 stop("Please choose a dNdS estimation method that is supported by this function.")
        
        if(!is.element(format,c("mase", "clustal", "phylip", "fasta" , "msf" )))
                stop("Please choose a format that is supported by seqinr::read.alignment.")
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- NULL
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        if(!file.exists(paste0("_calculation",f_sep))){
                
                dir.create("_calculation")
        }
        
        
        if(est.method == "Comeron"){
         
                if(is.null(subst_name)){
                        
                        file.out <- paste0("_calculation",f_sep,"gestimout")
                }
                
                if(!is.null(subst_name)){
                        
                        file.out <- paste0("_calculation",f_sep,subst_name,"_gestimout")    
                }
                
           # file in fasta required     
           if(format != "fasta")
                   stop("To use gestimator an alignment file in fasta format is required.")
           
            tryCatch(
            {    
                # To use gestimator a file in fasta format is required    

                # include this line instead of the following, to use internal gestimator
                gestimator(file = file, file_out=file.out) 
                #system(paste0("gestimator -i ",file," -o ",file.out))
                
                hit.table <-data.table::fread(file.out)
                data.table::setnames(hit.table, old=c("V1","V2","V3","V4","V5"), 
                                     new = c("query_id","subject_id","dN","dS","dNdS"))
                data.table::setkey(hit.table, query_id)
                
                if(!quiet){print("Substitutionrate successfully calculated by gestimator")}
                
                return(hit.table)
            },error = function(){ stop(paste0("Please check the correct path to ",est.method,
                                               "... the interface call did not work properly.") ) }
            
            )
        }
        
        if(est.method == "Li" ){
                
                        aln <- seqinr::read.alignment(file = file, format = format)
                
                        res_list <- seqinr::kaks(aln)
                        
                        res <- data.table::data.table(t(c(unlist(aln$nam),  res_list$ka, res_list$ks, (res_list$ka/res_list$ks))))
                        data.table::setnames(res, old = paste0("V",1:5), 
                                             new = c("query_id","subject_id", "dN", "dS","dNdS"))
                        data.table::setkey(res,query_id)
                 
                        if(!quiet){print("Substitutionrate successfully calculated using Li's method.")}
                        
                        return(res)
                        
        
        }
        
        if(is.element(est.method,kaks_calc_methods)){
                
                
                operating_sys <- Sys.info()[1]
                
                if (operating_sys == "Darwin"){ 
                        os <- "Mac"
                        calc <- "KaKs_Calculator"
                }
                
                if (operating_sys == "Linux"){
                        os <- "Linux"
                        calc <- "KaKs_Calculator"
                }
                
                if (operating_sys == "Windows"){ 
                        os <- "Windows"
                        calc <- "KaKs_Calculator.exe"
                }
                
                               
               # determine the file seperator of the current OS
               f_sep <- .Platform$file.sep
               fa2axt <- system.file(paste0("KaKs_Calculator1.2",f_sep,"parseFastaIntoAXT.pl"), package = "orthologr")
                              
               
               ### replace paths like this: path/to a folder/of_interest
               ### by: path/'to a folder'/of_interest
               file_name <- unlist(strsplit(file,f_sep))
               file_name <- file_name[length(file_name)]
               curr_wd <- unlist(strsplit(getwd(),f_sep))
               wdir <- grepl(" ",curr_wd)
               
               if(any(wdir)){
                       curr_wd[wdir] <- stringr::str_replace(string = curr_wd[wdir],replacement = paste0("'",curr_wd[wdir],"'"), pattern = curr_wd[wdir])
               }
               curr_wd <- paste0(curr_wd,collapse = f_sep)
               
               ###
               
               tryCatch(
                       
              {        # converting to input file for KaKs_Calculator
                       system(paste0("perl ",fa2axt ," ",file," ",curr_wd,f_sep,"_calculation",f_sep,file_name))
                       # running KaKs_Calculator inside the orthologr package
                       KaKs_Calculator <- system.file(paste0("KaKs_Calculator1.2",f_sep,"bin",f_sep,os,f_sep,calc), package = "orthologr")
               
                       if(is.null(kaks_calc.params))
                               system(paste0(KaKs_Calculator," -i ",paste0("_calculation",f_sep,file_name,".axt")," -o ",paste0("_calculation",f_sep,file_name,".axt.kaks"," -m ",est.method)))
               
                       if(!is.null(kaks_calc.params))
                               system(paste0(KaKs_Calculator," -i ",paste0("_calculation",f_sep,file_name,".axt")," -o ",paste0("_calculation",f_sep,file_name,".axt.kaks ",kaks_calc.params)))
               
              },error = function(){ stop(paste0("KaKs_Calculator 1.2 couln't run properly, please check your input files."))}
               
              )
              
              tryCatch(
                      {
                              kaks_tbl <- read.csv(paste0("_calculation",f_sep,file_name,".axt.kaks"),sep = "\t", header = TRUE)
                              kaks_tbl_res <- kaks_tbl[ , 1:5]
                              kaks_tbl_res <- data.frame(sapply(kaks_tbl_res[ , 1], function(x) unlist(strsplit(as.character(x),"-"))[1]),
                                                         sapply(kaks_tbl_res[ , 1], function(x) unlist(strsplit(as.character(x),"-"))[2]) ,
                                                         kaks_tbl_res[ , c(3:5,2)])
              
                              names(kaks_tbl_res) <- c("query_id","subject_id", "dN", "dS","dNdS","method")
                              kaks_tbl_res <- data.table::as.data.table(kaks_tbl_res)
                              data.table::setkey(kaks_tbl_res,query_id)
               
                              return(kaks_tbl_res)
                              
                      }, error = function(){stop(paste0("Something went wront with KaKs_Calculator .\n",
                                                         paste0("_calculation",f_sep,file_name,".axt.kaks"),
                                                         " could not be read properly."))}
              )
        }
        
}



#' @title This function computes the dNdS value of one given pairwise alignment.
#' @description This function takes a vector containing the
#' query amino acid sequence, subject amino acid sequence, query CDS sequence, and subject CDS sequence
#' and then runs the following pipieline:
#' 
#' 1) Multiple-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' @param complete_tbl a data.table object storing the query_id, subject_id, query_cds (sequence), 
#' subject_cds (sequence), query_aa (sequence), and subject_aa (sequence) of the organisms that shall be compared.
#' @param aa_aln_tool a character string specifying the multiple alignment tool that shall be used for pairwise protein alignments.
#' @param aa_aln_path a character string specifying the path to the corresponding multiple alignment tool.
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type} = "multiple" or \code{aa_aln_type} = "pairwise".
#' Default is \code{aa_aln_type} = "multiple".
#' @param aa_aln_params a character string specifying additional parameters that shall be passed to the multiple alignment system call.
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used for codon alignments. Default is \code{codon_aln_tool} = "pal2nal".
#' @param dnds_est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li", "YN", etc. See Details for all options.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @param comp_cores a numeric value specifying the number of cores that shall be used to perform parallel computations on a multicore machine.
#' @param clean_folders a boolean value spefiying whether all internall folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @author Hajk-Georg Drost and Sarah Scharfenberg
#' @details This function takes the amino acid and CDS sequences two orthologous genes
#' and writes the corresponding amino acid and CDS sequences as fasta file into
#' the internal folder environment. The resulting fasta files (two files) store the 
#' amino acid sequence of the query_id and subject_id (file one) and the CDS sequence of
#' the query_id and subject_id (file two). These fasta files are then used to pass through the following pipeline:
#' 
#' 1) Multiple-Alignment or Pairwise-Alignment of query amino acid sequence and subject amino acid sequence
#' 
#' 2) Codon-Alignment of the amino acid alignment returned by 1) and query CDS sequence + subject CDS sequence
#' 
#' 3) dNdS estimation of the codon alignment returned by 2)
#' 
#' @references http://www.r-bloggers.com/the-wonders-of-foreach/
#' @seealso \code{\link{multi_aln}}, \code{\link{substitutionrate}}, \code{\link{dNdS}}
#' @import foreach
#' @import data.table
compute_dnds <- function(complete_tbl,
                         aa_aln_type = "multiple", aa_aln_tool = "clustalw", aa_aln_path = NULL,
                         aa_aln_params = NULL, codon_aln_tool = "pal2nal",
                         dnds_est.method = "YN", quiet = FALSE, comp_cores = 1,
                         clean_folders = FALSE){
        
        if(comp_cores > parallel::detectCores())
                stop("You assigned more cores to the comp_cores argument than are availible on your machine.")
        
        
        # due to the discussion of no visible binding for global variable for
        # data.table objects see:
        # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
        query_id <- subject_id <- query_cds <- subject_cds <- query_aa <- subject_aa <- NULL
        
        multicore <- (comp_cores > 1)
        
        # determine the file seperator of the current OS
        f_sep <- .Platform$file.sep
        
        orthologs_names <- vector(mode = "list", length = 2)
        cds_seqs <- vector(mode = "list", length = 2)
        aa_seqs <- vector(mode = "list", length = 2)
        cds_session_fasta <- vector(mode = "character", length = 1)
        aa_session_fasta <- vector(mode = "character", length = 1)
        aa_aln_name <- vector(mode = "character", length = 1)
        
        if(!multicore)
                dNdS_values <- vector(mode = "list", length = nrow(complete_tbl))
        
        if(!file.exists(paste0("_alignment",f_sep))){
                
                dir.create("_alignment")
        }
        
        if(!file.exists(paste0("_alignment",f_sep,"orthologs",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"orthologs"))
        }
        
        if(!file.exists(paste0("_alignment",f_sep,"orthologs",f_sep,"CDS",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"orthologs",f_sep,"CDS"))
        }
        
        if(!file.exists(paste0("_alignment",f_sep,"orthologs",f_sep,"AA",f_sep))){
                
                dir.create(paste0("_alignment",f_sep,"orthologs",f_sep,"AA"))
        }
        
                
        if(!multicore){        
              for(i in 1:nrow(complete_tbl)){
                  dNdS_values[i] <- list((function(i) {
                
          ### Perform the sampling process in parallel
#         dNdS_values <- foreach::foreach(i = 1:nrow(complete_tbl),.combine="rbind") %dopar%{
                
                # storing the query gene id and subject gene id of the orthologous gene pair 
                orthologs_names <- list(complete_tbl[i , query_id],complete_tbl[i, subject_id])
        
                # storing the query CDS sequence and subject CDS sequence of the orthologous gene pair 
                cds_seqs <- list(complete_tbl[i, query_cds],complete_tbl[i, subject_cds])
        
                # storing the query amino acid sequence and subject amino acid sequence of the orthologous gene pair 
                aa_seqs <- list(complete_tbl[i, query_aa],complete_tbl[i, subject_aa])
        
        
                #aa_aln_name <- paste0("query_",i,"___","subject_",i)
                aa_aln_name <- paste0("q",i)

                # create cds fasta of orthologous gene pair having session name: 'aa_aln_name'
                cds_session_fasta <- paste0("_alignment",f_sep,"orthologs",f_sep,"CDS",f_sep,aa_aln_name,"_cds.fasta")
                seqinr::write.fasta(sequences = cds_seqs, names = orthologs_names, file.out = cds_session_fasta)
        
                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                aa_session_fasta <- paste0("_alignment",f_sep,"orthologs",f_sep,"AA",f_sep,aa_aln_name,"_aa.fasta")
                seqinr::write.fasta(sequences = aa_seqs, names = orthologs_names, file.out = aa_session_fasta)
        
                # which multi_aln tool should get the parameters
                #multi_aln_tool_params <- paste0(aa_aln_tool,".",params)
                          
                
#                 pairwise_aln <- Biostrings::pairwiseAlignment(aa_seqs[[1]],aa_seqs[[2]], type = "global")
#                 
#                 Biostrings::writePairwiseAlignments(pairwise_aln, block.width = 60)
#                 
                
                
              
                if(aa_aln_type == "multiple"){
                        
                        # align aa -> <aa_aln_tool>.aln
                        multi_aln(file = aa_session_fasta, 
                                  tool = aa_aln_tool, get_aln = FALSE, 
                                  multi_aln_name = aa_aln_name, 
                                  path = aa_aln_path, quiet = quiet)
        
                        aa_aln_session_name <- paste0(aa_aln_name,"_",aa_aln_tool,".aln")
                        
                        # align codon -> cds.aln
                        codon_aln(file_aln = paste0("_alignment",f_sep,"multi_aln",f_sep,aa_aln_session_name),
                                  file_nuc = cds_session_fasta, tool = codon_aln_tool,format = "fasta",
                                  codon_aln_name = aa_aln_name,
                                  get_aln = FALSE, quiet = quiet)
                }
                if(aa_aln_type == "pairwise"){
                        
                        # align aa -> <aa_aln_tool>.aln
                        pairwise_aln(file = aa_session_fasta, 
                                     tool = aa_aln_tool, 
                                     seq_type = "protein",
                                     get_aln = FALSE,
                                     pairwise_aln_name = aa_aln_name, 
                                     path = aa_aln_path, quiet = quiet)
                        
                        aa_aln_session_name <- paste0(aa_aln_name,"_",aa_aln_tool,"_AA.aln")
                        
                        # align codon -> cds.aln
                        codon_aln(file_aln = paste0("_alignment",f_sep,"pairwise_aln",f_sep,aa_aln_session_name),
                                  file_nuc = cds_session_fasta, tool = codon_aln_tool,format = "fasta",
                                  codon_aln_name = aa_aln_name,
                                  get_aln = FALSE, quiet = quiet)
                }
                
                
        
                codon_aln_session_name <- paste0(aa_aln_name,"_",codon_aln_tool,".aln")
                
                
                
                # compute kaks
                dNdS.table <- substitutionrate(file = paste0("_alignment",f_sep,"codon_aln",f_sep,codon_aln_session_name), 
                                               subst_name = aa_aln_name,
                                               est.method = dnds_est.method, quiet = quiet)
                
                return(dNdS.table)
        
                })(i)
              )
            }
        }


        if(multicore){
                ### Parallellizing the sampling process using the 'doMC' and 'parallel' package
                ### register all given cores for parallelization
                doMC::registerDoMC(comp_cores)
                
                ### Perform the sampling process in parallel
                dNdS_values <- foreach::foreach(i = 1:nrow(complete_tbl),.combine="rbind",
                                                .packages = c("seqinr","data.table"),
                                                .errorhandling = "stop") %dopar%{
                
                # storing the query gene id and subject gene id of the orthologous gene pair 
                orthologs_names <- list(complete_tbl[i , query_id],complete_tbl[i, subject_id])
                
                # storing the query CDS sequence and subject CDS sequence of the orthologous gene pair 
                cds_seqs <- list(complete_tbl[i, query_cds],complete_tbl[i, subject_cds])
                
                # storing the query amino acid sequence and subject amino acid sequence of the orthologous gene pair 
                aa_seqs <- list(complete_tbl[i, query_aa],complete_tbl[i, subject_aa])
                
                
                #aa_aln_name <- paste0("query_",i,"___","subject_",i)
                aa_aln_name <- paste0("q",i)
                
                # create cds fasta of orthologous gene pair having session name: 'aa_aln_name'
                cds_session_fasta <- paste0("_alignment",f_sep,"orthologs",f_sep,"CDS",f_sep,aa_aln_name,"_cds.fasta")
                seqinr::write.fasta(sequences = cds_seqs, names = orthologs_names, file.out = cds_session_fasta)
                
                # create aa fasta of orthologous gene pair having session name: 'aa_aln_name'
                aa_session_fasta <- paste0("_alignment",f_sep,"orthologs",f_sep,"AA",f_sep,aa_aln_name,"_aa.fasta")
                seqinr::write.fasta(sequences = aa_seqs, names = orthologs_names, file.out = aa_session_fasta)
                
                # which multi_aln tool should get the parameters
                #multi_aln_tool_params <- paste0(aa_aln_tool,".",params)
                
                #                 pairwise_aln <- Biostrings::pairwiseAlignment(aa_seqs[[1]],aa_seqs[[2]], type = "global")
                #                 
                #                 Biostrings::writePairwiseAlignments(pairwise_aln, block.width = 60)
                #       
                

                
                if(aa_aln_type=="multiple"){
                        
                        # align aa -> <aa_aln_tool>.aln
                        multi_aln(file = aa_session_fasta, 
                                  tool = aa_aln_tool, get_aln = FALSE, 
                                  multi_aln_name = aa_aln_name, 
                                  path = aa_aln_path, quiet = quiet)
                        
                        aa_aln_session_name <- paste0(aa_aln_name,"_",aa_aln_tool,".aln")
                        
                        # align codon -> cds.aln
                        codon_aln(file_aln = paste0("_alignment",f_sep,"multi_aln",f_sep,aa_aln_session_name),
                                  file_nuc = cds_session_fasta, tool = codon_aln_tool,format = "fasta",
                                  codon_aln_name = aa_aln_name,
                                  get_aln = FALSE, quiet = quiet)
                }
                
                if(aa_aln_type=="pairwise"){
                        
                        # align aa -> <aa_aln_tool>.aln
                        pairwise_aln(file = aa_session_fasta, 
                                     tool = aa_aln_tool, 
                                     seq_type = "protein",
                                     get_aln = FALSE,
                                     pairwise_aln_name = aa_aln_name, 
                                     path = aa_aln_path, quiet = quiet)
                        
                        aa_aln_session_name <- paste0(aa_aln_name,"_",aa_aln_tool,"_AA.aln")
                        
                        # align codon -> cds.aln
                        codon_aln(file_aln = paste0("_alignment",f_sep,"pairwise_aln",f_sep,aa_aln_session_name),
                                  file_nuc = cds_session_fasta, tool = codon_aln_tool,format = "fasta",
                                  codon_aln_name = aa_aln_name,
                                  get_aln = FALSE, quiet = quiet)
                }

                codon_aln_session_name <- paste0(aa_aln_name,"_",codon_aln_tool,".aln")
                
                # compute kaks
                dNdS.table <- substitutionrate(file = paste0("_alignment",f_sep,"codon_aln",f_sep,codon_aln_session_name), 
                                               subst_name = aa_aln_name,
                                               est.method = dnds_est.method, quiet = quiet)
                
                return(dNdS.table)
                
               }
                
        }


        if(!multicore)
                dNdS_tbl <- data.table::as.data.table(do.call(rbind,dNdS_values))

        if(multicore)
                dNdS_tbl <- data.table::as.data.table(dNdS_values)

        setkeyv(dNdS_tbl,c("query_id","subject_id"))
        

        if(clean_folders)
                clean_all_folders(c("_alignment", "_blast_db","_calculation"))

        # returning the dNdS table as data.table object
        return(dNdS_tbl)

}






