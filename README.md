orthologr
=========

### A framework combining common gene orthology detection methods and dNdS estimation models for pipeline processing 

The `orthologr` package provides interfaces between R and common bioinformatics tools
used to perform orthology inference, multiple sequence alignments, codon alignments, dNdS estimation,
CDS prediction, and divergence stratigraphy.

## Fast installation guide

```r
# install.packages("devtools")

# install the current version of orthologr on your system
library(devtools)
install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("orthologr", lib.loc = "C:/Program Files/R/R-3.1.1/library")


```


## Use Cases

Learn `orthologr` by use cases provided by the following tutorials: 

- [Orthology Inference using two or multiple species](https://github.com/HajkD/orthologr/blob/master/vignettes/orthology_inference.Rmd)
- [dNdS estimation for two or multiple species](https://github.com/HajkD/orthologr/blob/master/vignettes/dNdS_estimation.Rmd)
- [Phylotranscriptomics analyses based on sequence divergence: divergence stratigraphy](https://github.com/HajkD/orthologr/blob/master/vignettes/divergence_stratigraphy.Rmd)
- [Perform BLAST searches to find conserved regions in your set of subject organisms](https://github.com/HajkD/orthologr/blob/master/vignettes/blast.Rmd)
- [Sequence Alignments](https://github.com/HajkD/orthologr/blob/master/vignettes/sequence_alignments.Rmd)


## Installing prerequisite tools

Since many function have been implemented as interface functions between
the __R__ language and common bioinformatics tools, some of these tools must
be installed on your system to be able to use the corresponding functions in `orthologr`.

### Programming Languages

Most functions are optimized in terms of computational performance.
For this purpose, you need following programming languages on your system:

 - [__R__](http://www.cran.r-project.org) >= 3.1.1
 
 Please make sure you have the latest __R__ version installed on your system.
 
 To check which R version is currently installed on your system run the following command in __R__
 
```r
 # ckecking the current R version
 version$version.string
 
```
 
 - [__C++11__](http://isocpp.org/about)
 
 Please make also sure that you have C++11 installed on your system,
 this is important to use the features provided by [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html).
 
 - [__Perl__](https://www.perl.org) >= 5.12
 
 The Perl language is used to run [__Pal2Nal__](http://www.bork.embl.de/pal2nal/) and parts of
 [__KaKs_Calculator__](https://code.google.com/p/kaks-calculator/).
 
Run the following command in your Terminal to check whether you have [Perl](https://www.perl.org)
installed and what version is installed on your system
 
```shell
 perl -v
```
 
 - [__Python__](https://www.python.org) >= 2.7

Check the current [Python](https://www.python.org) version installed with

```shell
python --version

```



### Install additional R packages imported by orthologr

The `orthologr` package depends on a series of add ons (__R__ packages) that can
be downloaded as follows


```r
# install all CRAN packages on which orthologr depends on
install.packages(c("Rcpp","data.table","dplyr","doMC","foreach","ape"))

# install all Bioconductor packages orthologr depends on

# install Bioconductor base packages
source("http://bioconductor.org/biocLite.R")
biocLite()

# install package: Biostrings
biocLite("Biostrings")


```


### Orthology Inference Tools

The `orthologr` package provides interfaces to the following bioinformatics tools 
enabling orthology detection ([orthology inference](http://link.springer.com/protocol/10.1007%2F978-1-61779-582-4_9)):

The main function to perform orthology inference is `orthologs()` and internally this function relies
on the following interface functions implemented in `orthologr`.

 - [__BLAST__](http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/]) : Basic Local Alignment Search Tool finds regions of similarity between biological sequences and is also used as underlying paradigm of most fast orthology inference methods.
 
 There are several interface functions to BLAST+ implemented in `orthologr`
 
 - `blast()` : Interface function to BLAST+
 - `blast_best()` : Function to perform a BLAST+ best hit search
 - `blast_rec()` : Function to perform a BLAST+ reciprocal best hit (RBH) search
 - `set_blast()` : Function preparing the parameters and databases for subsequent BLAST+ searches
 - `advanced_blast()` : Advanced interface function to BLAST+
 - `advanced_makedb()` : Advanced interface function to makeblastdb
 
 - [__ProteinOrtho__](https://www.bioinf.uni-leipzig.de/Software/proteinortho/) : An orthology inference tool for large-scale analysis
 
 * `ProteinOrtho()`: Interface function to ProteinOrtho
 
 - [__OrthoMCL__](http://www.orthomcl.org/orthomcl/) : An Orthology Inference tool based on the OrthoMCL Algorithm detecting ortholog groups using all-versus all BLAST of all compared protein sequences
 
 * `OrthoMCL()`: Interface function to OrthoMCL
 
 
### Multiple Alignment Tools

The `orthologr` package also provides interfaces to the following Multiple Alignment Tools.
Nevertheless, non of them have to be installed if the corresponding interface functions
are not used.

 - [__ClustalW__](http://www.clustal.org/clustal2/) : Advanced multiple alignment tool of nucleic acid and protein sequences
 
 - [__T_Coffee__](http://www.tcoffee.org/Projects/tcoffee/) : A collection of tools for processing multiple sequence alignments
 of nucleic acids and proteins as well as their 3D structure
 
 - [__MUSCLE__](http://www.drive5.com/muscle/) : Fast and accurate multiple alignment tool of nucleic acid and protein sequences
 
 - [__ClustalO__](http://www.clustal.org/omega/) : Fast and scalable multiple alignment tool for nucleic acid and protein sequences that is also
 capable of performing HMM alignments
 
 - [__MAFFT__](http://mafft.cbrc.jp/alignment/software/) : A tool for multiple sequence alignment and phylogeny

In `orthologr` the function `multi_aln()` provides interfaces to all of these multiple alignment tools
as well as an pairwise alignment interface to the [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html) package performing a [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574).


### Codon Alignment Tools

The codon alignment tool [Pal2Nal](http://www.bork.embl.de/pal2nal/) is already integrated in the `orthologr` package
and doesn't need to be installed.

 - [__Pal2Nal__](http://www.bork.embl.de/pal2nal/)

You don't need to worry about downloading and installing __PAL2NAL__, it is already included in the `orthologr` package.
The corresponding function `codon_aln()` takes a protein alignment and the corresponding coding sequences and returns
a codon alignment by calling [__Pal2Nal__](http://www.bork.embl.de/pal2nal/) from inside of the `orthologr` package.

### dNdS Estimation Methods

dNdS estimation is a method to quantify the selection pressure acting on a specific protein sequence determined by pairwise comparisons of
amino acid substitutions between two protein sequences and their corresponding codon alignments.
Different models have been proposed to estimate this ratio quantifying selection pressure on proteins.
The `orthologr` package includes the most common dNdS estimation methods.

Starting with an codon alignment returned by `codon_aln()` the function `dNdS()` computes
the the dN, dS, and dNdS values of pairs of proteins.

Based on implementations provided by [gestimator](http://molpopgen.org/software/libsequence.html), [ape](http://www.cran.r-project.org/web/packages/ape/index.html), and [KaKs_Calculator](https://code.google.com/p/kaks-calculator/),
the following dNdS Estimation Methods are available in `orthologr`:

 - [Li](http://link.springer.com/article/10.1007/BF02407308#page-1) : Li's method (1993) -> provided by the [ape package](http://cran.r-project.org/web/packages/ape/index.html)

 - [Comeron](http://link.springer.com/article/10.1007/BF00173196) : Comeron's method (1995)

 - [NG](http://mbe.oxfordjournals.org/content/3/5/418.short) : Nei, M. and Gojobori, T. (1986)

 - [LWL](http://mbe.oxfordjournals.org/content/2/2/150.short) : Li, W.H., et al. (1985)

 - [MLWL](http://mbe.oxfordjournals.org/content/21/12/2290.short) (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)

 - [YN](http://mbe.oxfordjournals.org/content/17/1/32.short) : Yang, Z. and Nielsen, R. (2000)

 - [MYN](http://www.biomedcentral.com/1471-2148/6/44/) (Modified YN): Zhang, Z., et al. (2006)

For this purpose you need to have __KaKs_Calculator__ installed on your system and executable from your default `PATH`, e,g, `/usr/local/bin/`.



### CDS prediction pipeline

It is also possible to run __Divergence Stratigraphy__ or any other
orthology inference method starting with non-annotated genomes.

For this purpose `orthologr` provides a function `predict_cds()` that allows
you to perform gene prediction (CDS prediction) based on the gene prediction programs
[Augustus](http://bioinf.uni-greifswald.de/augustus/) and [GeneMark](http://exon.gatech.edu/GeneMark/)
and also allows you to translate the predicted genes to amino acids:  genes -> cds -> amino acids that
can subsequently be used for __Divergence Stratigraphy__ or __Orthology Inference__.

- [Augustus](http://bioinf.uni-greifswald.de/augustus/) : Prediction of regions within eukaryotic genomic sequences that can be denoted as "predicted genes"

* `augustus()`: Interface function for the gene prediction program augustus

- [GeneMark](http://exon.gatech.edu/GeneMark/): A family of gene prediction programs for all kingdoms of life

* `genemark()` : Interface function for the gene prediction program GeneMark



## Licenses

The `orthologr` package includes source code that has been published under following licenses:

### genevestigator

```
All files included in `orthologr` that were taken from genevestigator are 
also part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.


Modified by Sarah Scharfenberg and Hajk-Georg Drost (2014) to work 
in orthologr without using external libraries from libsequence.

All changes are also free under the terms of GNU General Public License
version 3 of the License, or any later version.

```

### KaKs_Calculator

In `orthologr` the file `parseFastaIntoAXT.pl` is stored in `/inst/KaKs_Calc_parser`.

```
The parseFastaIntoAXT.pl script is freely available under GNU GPL v3 
Licence and included in the KaKs_Calculator project that can be found at 
https://code.google.com/p/kaks-calculator/

```



