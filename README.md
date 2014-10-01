orthologr
=========

### A package combining common gene orthology detection methods and dNdS models for pipeline processing 

The `orthologr` package provides interfaces between R and common bioinformatics tools
used to perform orthology inference, multiple sequence alignments, codon alignments, dNdS computations,
and CDS predictions.

## Fast install guide

```r
# install.packages("devtools")

# install the current version of orthologr on your system
library(devtools)
install_github("HajkD/orthologr")

# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/orthologr")

# When working with Windows, first you need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/orthologr")

# and then call it from the library
library("orthologr", lib.loc = "C:/Program Files/R/R-3.1.1/library")

# when a github repo is set to 'private'
# you need to generate a authentification token and then use it via the argument: 'auth_token'

devtools::install_github("HajkD/orthologr", auth_token = "here your token")

```

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
 [KaKs_Calculator](https://code.google.com/p/kaks-calculator/).
 
Run the following command in your Terminal to check whether you have __Perl__
installed and what version is installed on your system
 
```shell
 perl -v
```
 
 - [__Python__](https://www.python.org) >= 2.8



### Install additional R packages imported by orthologr

The `orthologr` package depends on a series of add ons (__R__ packages) that can
be downloaded as follows


```r
# install all CRAN packages on which orthologr depends on
install.packages(c("Rcpp","data.table","dplyr","doMC","foreach","ape"))

# install all Bioconductor packages orthologr depends on
source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("Biostrings")


```


### Orthology Inference Tools

The `orthologr` package provides interfaces to the following bioinformatics tools 
enabling orthology detection (orthology inference):

 - [__BLAST__](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download])
 
 There are several interface functions to BLAST+ implemented in `orthologr`
 
 - `blast()` : Interface function to BLAST+
 - `blast_best()` : Function to perform a BLAST+ best hit search
 - `blast_rec()` : Function to perform a BLAST+ reciprocal best hit (RBH) search
 - `set_blast()` : Function preparing the parameters and databases for subsequent BLAST+ searches
 - `advanced_blast()` : Advanced interface function to BLAST+
 - `advanced_makedb()` : Advanced interface function to makeblastdb
 
 - [__ProteinOrtho__](https://www.bioinf.uni-leipzig.de/Software/proteinortho/)
 - [__OrthoMCL__](http://www.orthomcl.org/orthomcl/)
 - [__InParanoid__](http://inparanoid.sbc.su.se/cgi-bin/index.cgi)

### Multiple Alignment Tools

The `orthologr` package also provides interfaces to the following Multiple Alignment Tools.
Nevertheless, non of them have to be installed if the corresponding interface functions
are not used.

 - [__ClustalW__](http://www.clustal.org/clustal2/)
 - [__T_Coffee__](http://www.tcoffee.org/Projects/tcoffee/)
 - [__MUSCLE__](http://www.drive5.com/muscle/)
 - [__ClustalO__](http://www.clustal.org/omega/)
 - [__MAFFT__](http://mafft.cbrc.jp/alignment/software/)

In `orthologr` the function `multi_aln()` provides interfaces to all of these multiple alignment tools
as well as an pairwise alignment interface to the [Biostrings](http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html) package performing a Needleman-Wunsch algorithm.


### Codon Alignment Tools

The Codon Alignment Tools Pal2Nal is already integrated in the `orthologr` package
and doens't need to be installed.

 - [__Pal2Nal__](http://www.bork.embl.de/pal2nal/)

You don't need to worry about downloading and installing __PAL2NAL__, it is already included in the `orthologr` package.

### dNdS Estimation Methods

Based on implementations provided by [gestimator](http://molpopgen.org/software/libsequence.html), [ape](http://www.cran.r-project.org/web/packages/ape/index.html), and [KaKs_Calculator](https://code.google.com/p/kaks-calculator/),
the following dNdS Estimation Methods are available in `orthologr`:

 - "Li" : Li's method (1993) -> provided by the ape package

 - "Comeron" : Comeron's method (1995)

 - "NG": Nei, M. and Gojobori, T. (1986)

 - "LWL": Li, W.H., et al. (1985)

 - "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)

 - "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)

 - "YN": Yang, Z. and Nielsen, R. (2000)

 - "MYN" (Modified YN): Zhang, Z., et al. (2006)

You don't need to worry about downloading and installing __gestimator__ and __KaKs_Calculator__, they are already included in the `orthologr` package.



### CDS prediction pipeline







