#' @title Interface function for the gene prediction program augustus.
#' @description This function provides an interface to the gene prediction tool augustus.
#' 
#' @param file a character string specifying the path to the fasta file storing query genome.
#' @param species a character string specifying the species.
#' @param params a character string specifying additional arguments that are provided by augustus. Default is \code{params} = \code{NULL}.
#' A detailed set of additional parameters can be found here: \url{http://bioinf.uni-greifswald.de/augustus/binaries/README.TXT }.
#' @param file.out an optional character string spcifying the path where the output should
#' be saved. Default is \code{file.out} = \code{NULL}.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details The \code{augustus} function takes a string specifying the path to the genome file
#' of interest as first argument and the identifier for the species as second.
#' 
#' The input genome of interest must be provided as fasta file.
#'  
#' Genomes stored in fasta files can be downloaded from \url{http://ensemblgenomes.org/info/genomes}.
#'
#' 
#' species is one of the following identifiers
#' 
#' \tabular{lr}{
#' identifier          \tab species \cr
#' ------------------- \tab ---------------------\cr
#' human                                    \tab Homo sapiens\cr
#'         
#' fly                                      \tab Drosophila melanogaster\cr
#' 
#' arabidopsis                              \tab Arabidopsis thaliana\cr
#' 
#' brugia                                   \tab Brugia malayi\cr
#' 
#' aedes                                    \tab Aedes aegypti\cr
#' 
#' tribolium                                \tab Tribolium castaneum\cr
#' 
#' schistosoma                              \tab Schistosoma mansoni\cr
#' 
#' tetrahymena                              \tab Tetrahymena thermophila\cr
#' 
#' galdieria                                \tab Galdieria sulphuraria\cr
#' 
#' maize                                    \tab Zea mays\cr
#' 
#' toxoplasma                               \tab Toxoplasma gondii\cr
#' 
#' caenorhabditis                           \tab Caenorhabditis elegans\cr
#' 
#' (elegans)                                \tab Caenorhabditis elegans \cr
#' 
#' aspergillus_fumigatus                    \tab Aspergillus fumigatus\cr
#' 
#' aspergillus_nidulans                     \tab Aspergillus nidulans\cr
#' 
#' (anidulans)                              \tab Aspergillus nidulans\cr
#' 
#' aspergillus_oryzae                       \tab Aspergillus oryzae\cr
#' 
#' aspergillus_terreus                      \tab Aspergillus terreus\cr
#' 
#' botrytis_cinerea                         \tab Botrytis cinerea\cr
#' 
#' candida_albicans                         \tab Candida albicans\cr
#' 
#' candida_guilliermondii                   \tab Candida guilliermondii\cr
#' 
#' candida_tropicalis                       \tab Candida tropicalis\cr
#' 
#' chaetomium_globosum                      \tab Chaetomium globosum\cr
#' 
#' coccidioides_immitis                     \tab Coccidioides immitis\cr
#' 
#' coprinus                                 \tab Coprinus cinereus\cr
#' 
#' coprinus_cinereus                        \tab Coprinus cinereus\cr
#' 
#' coyote_tobacco                           \tab Nicotiana attenuata\cr
#' 
#' cryptococcus_neoformans_gattii           \tab Cryptococcus neoformans gattii\cr
#' 
#' cryptococcus_neoformans_neoformans_B     \tab Cryptococcus neoformans neoformans\cr
#' 
#' cryptococcus_neoformans_neoformans_JEC21 \tab Cryptococcus neoformans neoformans\cr
#' 
#' (cryptococcus)                           \tab Cryptococcus neoformans\cr
#' 
#' debaryomyces_hansenii                    \tab Debaryomyces hansenii\cr
#' 
#' encephalitozoon_cuniculi_GB              \tab Encephalitozoon cuniculi\cr
#' 
#' eremothecium_gossypii                    \tab Eremothecium gossypii\cr
#' 
#' fusarium_graminearum                     \tab Fusarium graminearum\cr
#' 
#' (fusarium)                               \tab Fusarium graminearium\cr
#' 
#' histoplasma_capsulatum                   \tab Histoplasma capsulatum\cr
#' 
#' (histoplasma)                            \tab Histoplasma capsulatum\cr
#' 
#' kluyveromyces_lactis                     \tab Kluyveromyces lactis\cr
#' 
#' laccaria_bicolor                         \tab Laccaria bicolor\cr
#' 
#' lamprey                                  \tab Petromyzon marinus\cr
#' 
#' leishmania_tarentolae                    \tab Leishmania tarentolae\cr
#' 
#' lodderomyces_elongisporus                \tab Lodderomyces elongisporus\cr
#' 
#' magnaporthe_grisea                       \tab Magnaporthe grisea\cr
#' 
#' neurospora_crassa                        \tab Neurospora crassa\cr
#' 
#' (neurospora)                             \tab Neurospora crassa\cr
#' 
#' phanerochaete_chrysosporium              \tab Phanerochaete chrysosporium\cr
#' 
#' (pchrysosporium)                         \tab Phanerochaete chrysosporium\cr
#' 
#' pichia_stipitis                          \tab Pichia stipitis\cr
#' 
#' rhizopus_oryzae                          \tab Rhizopus oryzae\cr
#' 
#' saccharomyces_cerevisiae_S288C           \tab Saccharomyces cerevisiae\cr
#' 
#' saccharomyces_cerevisiae_rm11-1a_1       \tab Saccharomyces cerevisiae\cr
#' 
#' (saccharomyces)                          \tab Saccharomyces cerevisiae\cr
#' 
#' schizosaccharomyces_pombe                \tab Schizosaccharomyces pombe\cr
#' 
#' thermoanaerobacter_tengcongensis         \tab Thermoanaerobacter tengcongensis\cr
#' 
#' trichinella                              \tab Trichinella spiralis\cr
#' 
#' ustilago_maydis                          \tab Ustilago maydis\cr
#' 
#' (ustilago)                               \tab Ustilago maydis\cr
#' 
#' yarrowia_lipolytica                      \tab Yarrowia lipolytica\cr
#' 
#' nasonia                                  \tab Nasonia vitripennis\cr
#' 
#' tomato                                   \tab Solanum lycopersicum\cr
#' 
#' chlamydomonas                            \tab Chlamydomonas reinhardtii\cr
#' 
#' amphimedon                               \tab Amphimedon queenslandica\cr
#' 
#' pneumocystis                             \tab Pneumocystis jirovecii\cr
#' 
#' wheat                                    \tab Triticum aestivum\cr
#' 
#' chicken                                  \tab Gallus gallus\cr
#' 
#' zebrafish                                \tab Danio rerio\cr
#' 
#' E_coli_K12                               \tab Escherichia coli\cr
#' 
#' s_aureus                                 \tab Staphylococcus aureus
#'
#'}
#' @references 
#' 
#' Oliver Keller, Martin Kollmar, Mario Stanke, Stephan Waack (2011)
#' A novel hybrid gene prediction method employing protein multiple sequence alignments
#'Bioinformatics, doi: 10.1093/bioinformatics/btr010
#'
#'Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008)
#'Using native and syntenically mapped cDNA alignments to improve de novo gene finding
#'Bioinformatics, doi: 10.1093/bioinformatics/btn013
#'
#'Mario Stanke, Ana Tzvetkova, Burkhard Morgenstern (2006)
#'"AUGUSTUS at EGASP: using EST, protein and genomic alignments for improved gene prediction in the human genome"
#'BMC Genome Biology, 7(Suppl 1):S11.
#'
#'M. Stanke , O. Schoeffmann , B. Morgenstern, S. Waack (2006)
#'Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources 
#'BMC Bioinformatics 7, 62.
#'
#'Mario Stanke and Burkhard Morgenstern (2005)
#'"AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints",
#'Nucleic Acids Research, 33, W465-W467
#'
#'Mario Stanke, Rasmus Steinkamp, Stephan Waack and Burkhard Morgenstern (2004) 
#'"AUGUSTUS: a web server for gene finding in eukaryotes" 
#'Nucleic Acids Research, Vol. 32, W309-W312
#'
#'Mario Stanke and Stephan Waack (2003)
#'Gene Prediction with a Hidden-Markov Model and a new Intron Submodel. 
#'Bioinformatics, Vol. 19, Suppl. 2, pages ii215-ii225
#'Mario Stanke (2003)
#'
#'Gene Prediction with a Hidden-Markov Model. 
#'Ph.D. thesis, Universitaet Goettingen
#'
#' augustus: \url{http://bioinf.uni-greifswald.de/augustus/}
#' 
#' download augustus: \url{http://stubber.math-inf.uni-greifswald.de/bioinf/augustus/binaries/}
#' 
#'  @examples \dontrun{
#'  
#'  
#' # run example with augustus
#' augustus( file=system.file('seqs/example.fa', package = 'orthologr'), 
#' species="human", params = "--UTR=on") 
#' 
#' }
#'
#' @return augustus returns...
#' @export
augustus <- function(file, species, params = NULL, file.out = NULL){
        
        tryCatch(
                {
                        # augustus [parameters] --species=SPECIES queryfilename
                        if(is.null(file.out)){
                                
                                if(!is.null(params))
                                        system(paste0("augustus ",params," --species=",species," ",file))
                                
                                if(is.null(params))
                                        system(paste0("augustus --species=",species," ",file))
                                
                        } else {
                                
                                if(!is.null(params))
                                        system(paste0("augustus ",params," --species=",species," ",file," >>",file.out))
                                
                                if(is.null(params))
                                        system(paste0("augustus --species=",species," ",file," >>",file.out))
                                
                                
                        }  
        
                }, error = function() stop("The augustus interface did not work properly, please make sure you have augustus installed properly and
                                           that it can be called from you default executable PATH.")
        )
}

is.provided_species <- function(species){
        return(is.element(species, 
            c("fly", "arabidopsis", "brugia", "aedes", "tribolium", "schistosoma", 
              "tetrahymena", "galdieria", "maize", "toxoplasma", "caenorhabditis", 
              "(elegans)", "aspergillus_fumigatus", "aspergillus_nidulans",
              "(anidulans)", "aspergillus_oryzae", "aspergillus_terreus", 
              "botrytis_cinerea", "candida_albicans", "candida_guilliermondii", 
              "candida_tropicalis", "chaetomium_globosum", "coccidioides_immitis", 
              "coprinus", "coprinus_cinereus", "coyote_tobacco", 
              "cryptococcus_neoformans_gattii", "cryptococcus_neoformans_neoformans_B", 
              "cryptococcus_neoformans_neoformans_JEC21", "(cryptococcus)", 
              "debaryomyces_hansenii", "encephalitozoon_cuniculi_GB", 
              "eremothecium_gossypii", "fusarium_graminearum", "(fusarium)", 
              "histoplasma_capsulatum", "(histoplasma)", "kluyveromyces_lactis", 
              "laccaria_bicolor", "lamprey", "leishmania_tarentolae", 
              "lodderomyces_elongisporus", "magnaporthe_grisea", "neurospora_crassa", 
              "(neurospora)", "phanerochaete_chrysosporium", "(pchrysosporium)",
              "pichia_stipitis", "rhizopus_oryzae", "saccharomyces_cerevisiae_S288C",
              "saccharomyces_cerevisiae_rm11-1a_1", "(saccharomyces)",
              "schizosaccharomyces_pombe", "thermoanaerobacter_tengcongensis", 
              "trichinella", "ustilago_maydis", "(ustilago)", "yarrowia_lipolytica", 
              "nasonia", "tomato", "chlamydomonas", "amphimedon", "pneumocystis",
              "wheat", "chicken", "zebrafish", "E_coli_K12", "s_aureus")))
}
