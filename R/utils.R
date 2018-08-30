test <- function(x){ print(paste0("Test ",x," passed.","\n"))}


is.dnaSequence <- function(seq){
        
        return((length(grep("[^(A|C|G|T|a|c|g|t)]", seq)) == 0))
}


# This functions are internal functions to proof if the used tool is provided. 

is.ortho_detection_method <- function(method = NULL){
        
        provided <- c("BH","RBH","PO","DELTA")
        
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
        kaks_calc_methods <- c("MA","MS","NG","LWL","LPB","MLWL","YN","MYN","GY","GMYN","ALL","kaks_calc")
        
        provided <- c("Comeron","Li",kaks_calc_methods)
        
        if (is.null(method)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(method, provided))
}

is.alignment_search_tool <- function(tool = NULL){
        
        provided <-  c("ggsearch","ssearch")
        
        if (is.null(tool)){
                print("The followig methods are provided: ")
                print(provided)
                return(FALSE)
        }
        
        return(is.element(tool, provided))
}

#' @title Translate CDS file to Amino Acids file
#' @description This function takes a \code{data.table} object as input and
#' translates the cds sequences stored as \code{seqs} column into the corresponding amino acid
#' sequence.
#' @param dt a \code{data.table} object storing cds sequences in a column named \code{seqs}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @import data.table
cds2aa <- function(dt, delete_corrupt_cds = FALSE){
        
        if (!is.data.table(dt))
                stop("Your CDS file was not corretly transformed into a data.table object.")
        
        # define visible bindings for global variables
        seqs <- aa <- geneids <- NULL
        
        # When using data.tables within packaes, always make sure
        # 'data.table' is included in the DESCRIPTION file as 'Imports' AND
        # in the NAMESPACE file with 'imports(data.table)' -> @import when using roxygen2
        # http://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package
        # https://github.com/hadley/dplyr/issues/548
        
        # omit empty sequences
        dt <- dt[ , .SD[sapply(seqs, function(x) {
                return(!(is.na(x) || x == ""))
        })]]
        
        if (delete_corrupt_cds) {
                # omit sequences that are not multiples of 3
                dt <- dt[ , .SD[sapply(seqs, function(x) {
                        return(nchar(x) %% 3 == 0)
                })]]
        }
        
        # omit sequences consisting of others than ACGT
        dt <- dt[ , .SD[sapply(seqs, is.dnaSequence)]]
        
        # translate cds to protein sequences
        tryCatch({
                dt[ , aa := transl(seqs), by = geneids]
                
        }, error = function(e) {
                stop(
                        "The input coding sequences could not be translated properly to amino acid sequences.",
                        "\n",
                        " Please check whether ",
                        file,
                        " stores valid coding sequences."
                )
        })
        
        return(dt)
}
        

