#' @title Compute a Codon Alignment
#' @description This function takes a protein alignment and the corresponding coding sequences 
#' as fasta files and computes the corresponding codon alignment based on the models specified in the PAL2NAL program.
#' @param file_aln a character string specifying the path to the file storing the protein alignment in CLUSTAL or FASTA format.
#' @param file_nuc a character string specifying the path to the file storing the coding sequences in multiple FASTA format.
#' @param format a character string specifying the file format used to store the codon alignment, e.g. "fasta", "clustal".
#' @param tool a character string specifying the program that should be used e.g. "pal2nal". 
#' @param codon_aln_name a character string specifying the name of the stored alignment file. 
#' @param params a character string specifying additional parameters that shall be passed to PAL2NAL. For example
#' \code{params} = \code{"-codontable 2"} would specify a codon table: \emph{vertebrate mitochondrial code}
#' instead if the \emph{universal code} (default). The default value is \code{params} = \code{NULL}.
#' Default is \code{codon_aln_name} = \code{NULL} denoting a default name: 'tool_name.aln' .
#' @param get_aln a logical value indicating whether the produced alignment should be returned.
#' @param quiet a logical value specifying whether a successful interface call shall be printed out.
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function provides an interface between the R language and common codon alignment tools such as "PAL2NAL".
#' Codon alignments can be used to quantify the evolutionary pressure acting on protein sequences based on dNdS estimation. 
#' 
#' The PAL2NAL program is used to convert a multiple sequence alignment or pairwise alignment of proteins and the corresponding DNA (or mRNA) sequences into 
#' a codon-based DNA alignment. The output codon-based DNA alignment can then be used for \code{\link{dNdS}} estimation.
#' 
#' The popular codon alignment tool PAL2NAL is included in this package and is being called when choosing \code{tool} = \code{"pal2nal"}.
#' @note
#' 
#' The PAL2NAL program automatically checks
#' \itemize{
#' \item If you use the same IDs are used in the protein alignment file and DNA (or mRNA) file ->
#' in this case the sequences don't have to be in the same order
#' 
#' \item If you don't use the same IDs are used in the protein alignment file and DNA (or mRNA) file
#' -> you have to rearrange the IDs or sequences so that sequences in the protein alignment file and DNA (or mRNA) file
#' have the same order
#' 
#' }
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' # performing a codon alignment using PAL2NAL
#' codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#'                        file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'), 
#'                        format   = "clustal", 
#'                        tool     = "pal2nal", 
#'                        get_aln  = TRUE)
#'                        
#'                        
#' # running PAL2NAL with additional parameters (e.g. a different codon table)
#' codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#'                        file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'), 
#'                        format   = "clustal", 
#'                        tool     = "pal2nal", 
#'                        get_aln  = TRUE,
#'                        params   = "-codontable 2")
#'                         
#' }
#'                        
#' @return In case \code{get_aln} = \code{TRUE} an seqinr alignment object is returned.
#' @references 
#' 
#' Mikita Suyama, David Torrents, and Peer Bork (2006)
#' PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Res. 34, W609-W612.
#' 
#' http://www.bork.embl.de/pal2nal/
#' 
#' http://www.genome.med.kyoto-u.ac.jp/cgi-bin/suyama/pal2nal/index.cgi
#' 
#' http://abacus.gene.ucl.ac.uk/software/paml.html
#' 
#' @seealso \code{\link{pairwise_aln}}, \code{\link{multi_aln}}, \code{\link{substitutionrate}},
#'  \code{\link{dNdS}}, \code{\link{divergence_stratigraphy}}
#' @export
codon_aln <- function(file_aln, 
                      file_nuc, 
                      format         = "clustal", 
                      tool           = "pal2nal", 
                      params         = NULL, 
                      codon_aln_name = NULL, 
                      get_aln        = FALSE, 
                      quiet          = FALSE){
        
        
        # the Pal2Nal program is stored as executable within
        # the R package environment: 'exec' folder
        # this is not an elegant version, only motivated by this discussion:
        # http://stackoverflow.com/questions/13463103/inst-and-extdata-folders-in-r-packaging
        # is there a better choice?
        path <-
                system.file(
                        file.path("pal2nal", "pal2nal.v14", "pal2nal.pl"),
                        package = "orthologr",
                        mustWork = TRUE
                )
        #path <- "/exec/pal2nal.v14/"
        
        if (!is.element(tool, c("pal2nal")))
                stop("Please choose a tool that is supported by this function.")
        
        if (!is.element(format, c("clustal", "fasta")))
                stop("Please choose a format that is supported by this function.")
        
        if (!file.exists(file.path(tempdir(), "_alignment"))) {
                dir.create(file.path(tempdir(), "_alignment"))
        }
        
        
        if (!file.exists(file.path(tempdir(), "_alignment", "codon_aln"))) {
                dir.create(file.path(tempdir(), "_alignment", "codon_aln"))
        }
        
        if (is.null(codon_aln_name))
                file.out <-
                file.path(tempdir(), "_alignment", "codon_aln", paste0(tool, ".aln"))
        
        if (!is.null(codon_aln_name))
                file.out <-
                file.path(tempdir(),
                          "_alignment",
                          "codon_aln",
                          paste0(codon_aln_name, "_", tool, ".aln"))
        
        
        # test whether the connection to pal2nal works
        
        if (tool == "pal2nal") {
                
                operating_sys <- Sys.info()[1]
                
                if (operating_sys == "Windows") {
                        tryCatch({
                                if (is.null(params))
                                        shell(
                                                paste0(
                                                        "perl ",
                                                        path,
                                                        " ",
                                                        file_aln,
                                                        " ",
                                                        file_nuc,
                                                        " -output ",
                                                        format,
                                                        " >",
                                                        file.out
                                                )
                                        )
                                
                                
                                if (!is.null(params))
                                        shell(
                                                paste0(
                                                        "perl ",
                                                        path,
                                                        " ",
                                                        file_aln,
                                                        " ",
                                                        file_nuc,
                                                        " -output ",
                                                        format,
                                                        " >",
                                                        file.out,
                                                        " ",
                                                        params
                                                )
                                        )
                                
                                
                        }, error = function(e) {
                                stop(
                                        "Please check the correct path to ",
                                        tool,
                                        "... the interface call did not work properly."
                                )
                        })
                } else {
                
                tryCatch({
                        if (is.null(params))
                                system(
                                        paste0(
                                                "perl ",
                                                path,
                                                " ",
                                                file_aln,
                                                " ",
                                                file_nuc,
                                                " -output ",
                                                format,
                                                " >",
                                                file.out
                                        )
                                )
                        
                        
                        if (!is.null(params))
                                system(
                                        paste0(
                                                "perl ",
                                                path,
                                                " ",
                                                file_aln,
                                                " ",
                                                file_nuc,
                                                " -output ",
                                                format,
                                                " >",
                                                file.out,
                                                " ",
                                                params
                                        )
                                )
                        
                        
                }, error = function(e) {
                        stop(
                                "Please check the correct path to ",
                                tool,
                                "... the interface call did not work properly."
                        )
                })
            }
        }
        
        if (!quiet) {
                print(paste0("Codon Alignment successfully written in ", file.out, "."))
        }
        
        
        if (get_aln) {
                tryCatch({
                        dna_aln <-
                                seqinr::read.alignment(file = file.out, format = "clustal")
                        return(dna_aln)
                        
                }, error = function(e) {
                        stop(
                                "Something went wront with Pal2Nal.pl .",
                                "\n",
                                file.out,
                                " could not be read properly."
                        )
                })
        }
        
}











