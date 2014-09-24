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

```

## Installing prerequisite tools

Since many function have been implemented as interface functions between
the __R__ language and common bioinformatics tools, some of these tools must
be installed on your system to be able to use the corresponding functions in `orthologr`.

### Programming Languages

Most functions are optimized in terms of computational performance.
For this purpose, you need following programming languages on your system:

#### - C++11
#### - Perl >= 5.12
#### - Python >= 2.8
#### - R >= 3.1.1


### Orthology Inference Tools

The `orthologr` package provides interfaces to the following bioinformatics tools 
enabling orthology detection (orthology inference):

#### - BLAST
#### - ProteinOrtho
#### - OrthoMCL
#### - InParanoid

### Multiple Alignment Tools

The `orthologr` package also provides interfaces to the following Multiple Alignment Tools.
Nevertheless, non of them have to be installed if the corresponding interface functions
are not used.

#### - ClustalW
#### - T_Coffee
#### - MUSCLE
#### - ClustalO
#### - MAFFT


### Codon Alignment Tools

The Codon Alignment Tools Pal2Nal is already integrated in the `orthologr` package
and doens't need to be installed.

#### - Pal2Nal

### dNdS Estimation Methods

Based on implementations provided by gestimator, ape, and KaKs_Calculator,
the following dNdS Estimation Methods are available in `orthologr`:

##### - "Li" : Li's method (1993) -> provided by the ape package

##### - "Comeron" : Comeron's method (1995)

##### - "NG": Nei, M. and Gojobori, T. (1986)

##### - "LWL": Li, W.H., et al. (1985)

##### - "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)

##### - "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)

##### - "YN": Yang, Z. and Nielsen, R. (2000)

##### - "MYN" (Modified YN): Zhang, Z., et al. (2006)





