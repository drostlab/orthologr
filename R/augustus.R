#' @title Wrapper for the gene prediction program augustus.
#' @description This function is a wrapper of the gene prediction tool augustus that can 
#' be downloaded here: http://stubber.math-inf.uni-greifswald.de/bioinf/augustus/binaries/ .
#' @param file a character string specifying the path to the fasta file storing query genome.
#' @param species a character string specifying the species.
#' @param params a character string specifying additional arguments that are provided by augustus.
#' @param file.out an optional character string spcifying the path where the output should
#' be saved.
#' @author Sarah Scharfenberg and Hajk-Georg Drost
#' @details The \code{augustus} function takes a string specifying the path to the genome file
#' of interest as first argument and the identifier for the species as second.
#' 
#' The genome ir required in fasta format. 
#' Genomes stored in fasta files can be downloaded from http://ensemblgenomes.org/info/genomes.
#'
#' Additional parameters can be added as a string and can be found here
#' http://bioinf.uni-greifswald.de/augustus/binaries/README.TXT .
#' 
#' species is one of the following identifiers
#' 
#' identifier                               | species
#' -----------------------------------------|----------------------
#'         human                                    | Homo sapiens
#' fly                                      | Drosophila melanogaster
#' arabidopsis                              | Arabidopsis thaliana
#' brugia                                   | Brugia malayi
#' aedes                                    | Aedes aegypti
#' tribolium                                | Tribolium castaneum
#' schistosoma                              | Schistosoma mansoni
#' tetrahymena                              | Tetrahymena thermophila
#' galdieria                                | Galdieria sulphuraria
#' maize                                    | Zea mays
#' toxoplasma                               | Toxoplasma gondii
#' caenorhabditis                           | Caenorhabditis elegans
#' (elegans)                                | Caenorhabditis elegans 
#' aspergillus_fumigatus                    | Aspergillus fumigatus
#' aspergillus_nidulans                     | Aspergillus nidulans
#' (anidulans)                              | Aspergillus nidulans
#' aspergillus_oryzae                       | Aspergillus oryzae
#' aspergillus_terreus                      | Aspergillus terreus
#' botrytis_cinerea                         | Botrytis cinerea
#' candida_albicans                         | Candida albicans
#' candida_guilliermondii                   | Candida guilliermondii
#' candida_tropicalis                       | Candida tropicalis
#' chaetomium_globosum                      | Chaetomium globosum
#' coccidioides_immitis                     | Coccidioides immitis
#' coprinus                                 | Coprinus cinereus
#' coprinus_cinereus                        | Coprinus cinereus
#' coyote_tobacco                           | Nicotiana attenuata
#' cryptococcus_neoformans_gattii           | Cryptococcus neoformans gattii
#' cryptococcus_neoformans_neoformans_B     | Cryptococcus neoformans neoformans
#' cryptococcus_neoformans_neoformans_JEC21 | Cryptococcus neoformans neoformans
#' (cryptococcus)                           | Cryptococcus neoformans
#' debaryomyces_hansenii                    | Debaryomyces hansenii
#' encephalitozoon_cuniculi_GB              | Encephalitozoon cuniculi
#' eremothecium_gossypii                    | Eremothecium gossypii
#' fusarium_graminearum                     | Fusarium graminearum
#' (fusarium)                               | Fusarium graminearium
#' histoplasma_capsulatum                   | Histoplasma capsulatum
#' (histoplasma)                            | Histoplasma capsulatum
#' kluyveromyces_lactis                     | Kluyveromyces lactis
#' laccaria_bicolor                         | Laccaria bicolor
#' lamprey                                  | Petromyzon marinus
#' leishmania_tarentolae                    | Leishmania tarentolae
#' lodderomyces_elongisporus                | Lodderomyces elongisporus
#' magnaporthe_grisea                       | Magnaporthe grisea
#' neurospora_crassa                        | Neurospora crassa
#' (neurospora)                             | Neurospora crassa
#' phanerochaete_chrysosporium              | Phanerochaete chrysosporium
#' (pchrysosporium)                         | Phanerochaete chrysosporium
#' pichia_stipitis                          | Pichia stipitis
#' rhizopus_oryzae                          | Rhizopus oryzae
#' saccharomyces_cerevisiae_S288C           | Saccharomyces cerevisiae
#' saccharomyces_cerevisiae_rm11-1a_1       | Saccharomyces cerevisiae
#' (saccharomyces)                          | Saccharomyces cerevisiae
#' schizosaccharomyces_pombe                | Schizosaccharomyces pombe
#' thermoanaerobacter_tengcongensis         | Thermoanaerobacter tengcongensis
#' trichinella                              | Trichinella spiralis
#' ustilago_maydis                          | Ustilago maydis
#' (ustilago)                               | Ustilago maydis
#' yarrowia_lipolytica                      | Yarrowia lipolytica
#' nasonia                                  | Nasonia vitripennis
#' tomato                                   | Solanum lycopersicum
#' chlamydomonas                            | Chlamydomonas reinhardtii
#' amphimedon                               | Amphimedon queenslandica
#' pneumocystis                             | Pneumocystis jirovecii
#' wheat                                    | Triticum aestivum
#' chicken                                  | Gallus gallus
#' zebrafish                                | Danio rerio
#' E_coli_K12                               | Escherichia coli
#' s_aureus                                 | Staphylococcus aureus
#'
#'  @examples \dontrun{
#' # run example with augustus
#' augustus( file="/home/sarah/Programs/augustus-3.0.3/examples/example.fa", 
#' species="human", params = "--UTR=on") 
#' }
#'
#' @return I dont know
#' @export
augustus <- function(file, species, file.out=NULL, params=""){
                
        # augustus [parameters] --species=SPECIES queryfilename
        if(is.null(file.out)){
                system(paste0("augustus ",params," --species=",species," ",file))
        }
        else {
                system(paste0("augustus ",params," --species=",species," ",file," >>",file.out))
        }  
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
