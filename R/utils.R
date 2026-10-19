# Suppress R CMD CHECK NOTEs for data.table's `.()` alias and other NSE symbols
utils::globalVariables(".")

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
cds2aa <- function(dt, delete_corrupt_cds = TRUE){
        
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
        



# internal: make sure `path` is on PATH for the lifetime of the R session,
# then return "" so this can be spliced into paste0() calls that build shell
# commands without changing their structure.
#
# Earlier versions of this helper returned a shell-syntax prefix such as
# 'export PATH=$PATH:<path>; ' to prepend to the command string. That is
# POSIX shell syntax and breaks on Windows, where system() may dispatch to
# cmd.exe (no 'export', and 'set VAR=val && cmd' doesn't see the new VAR
# because cmd.exe expands %VAR% before 'set' runs) or to an sh.exe picked up
# from Rtools/Git on PATH. Setting PATH in the R process itself works
# identically everywhere, since child processes always inherit R's
# environment.
path_prefix <- function(path) {
        if (is.null(path))
                return("")
        
        path <- normalizePath(path, mustWork = FALSE)
        current <- Sys.getenv("PATH")
        current_entries <- strsplit(current, .Platform$path.sep, fixed = TRUE)[[1]]
        
        if (!(path %in% current_entries))
                Sys.setenv(PATH = paste(path, current, sep = .Platform$path.sep))
        
        ""
}

# internal: run a shell command that relies on '>' or '2>' shell redirection
# (e.g. to capture a tool's stdout into an alignment output file). Plain
# system() does not reliably route through cmd.exe on every Windows setup
# (e.g. when an sh.exe from Rtools/Git is picked up on PATH), which silently
# breaks redirection: the literal '>file' text ends up passed to the program
# as just another argument instead of being interpreted by a shell. shell()
# explicitly always invokes the Windows command interpreter.
system_with_redirect <- function(command) {
        if (.Platform$OS.type == "windows") {
                suppressWarnings(shell(command, mustWork = FALSE))
        } else {
                system(command)
        }
}

# internal: run a shell command and return its exit status.
# If quiet = TRUE, stderr is captured instead of printed (e.g. DIAMOND's harmless
# "Failed to delete temporary file" warnings) and only shown if the command fails.
run_system_quiet <- function(command, quiet = TRUE) {
        if (!quiet)
                return(system(command))
        
        err_file <- tempfile("stderr_", fileext = ".txt")
        on.exit(unlink(err_file), add = TRUE)
        
        redirected <- if (.Platform$OS.type == "windows") {
                paste0(command, " 2>", shQuote(err_file))
        } else {
                # group the command in case it is itself a compound shell command
                paste0("{ ", command, "\n} 2>", shQuote(err_file))
        }
        
        # plain system() does not reliably route through cmd.exe on every
        # Windows setup (e.g. when an sh.exe from Rtools/Git is picked up),
        # which silently breaks '2>' redirection: the literal text ends up
        # passed to the program as just another argument. shell() explicitly
        # always invokes the Windows command interpreter, so redirection works.
        status <- if (.Platform$OS.type == "windows") {
                suppressWarnings(shell(redirected, mustWork = FALSE))
        } else {
                system(redirected)
        }
        
        if (!identical(status, 0L) && file.exists(err_file)) {
                err <- readLines(err_file, warn = FALSE)
                if (length(err) > 0)
                        message(paste(err, collapse = "\n"))
        }
        status
}
