test <- function(x){ print(paste0("Test ",x," passed.","\n"))}


is.dnaSequence <- function(seq){
        
        return((length(grep("[^(A|C|G|T|a|c|g|t)]", seq)) == 0))
}


# This functions are internal functions to proof if the used tool is provided. 

is.ortho_detection_method <- function(method = NULL){
        
        provided <- c("BH","RBH","PO","OrthoMCL","DELTA","GGSEARCH","SSEARCH")
        
        if(is.null(method)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(method, provided))
}


is.multiple_aln_tool <- function(tool = NULL){
        
        provided <- c("clustalw", "clustalo","muscle", "t_coffee", "mafft")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.pairwise_aln_tool <- function(tool = NULL){
        
        provided <- c("NW")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.codon_aln_tool <- function(tool = NULL){
        
        provided <-  c("pal2nal")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

is.dnds_est_method <- function(method = NULL){
        
        # dNdS estimation methods provided by the KaKs_Calculator 1.2 program
        kaks_calc_methods <- c("MA","MS","NG","LWL","LPB","MLWL","YN","MYN","GY","kaks_calc")
        
        provided <- c("Comeron","Li",kaks_calc_methods)
        
        if(is.null(method)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(method, provided))
}

is.alignment_search_tool <- function(tool = NULL){
        
        provided <-  c("ggsearch","ssearch")
        
        if(is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}


#         
#         is. ... <- function(tool=NULL){
#                 provided <-  c( ... )
#                 if(is.null(tool)){
#                         print("The followig methods are provided: ")
#                         print(provided)
#                         return(FALSE)
#                 }
#                 return(is.element(tool, provided))
#         }