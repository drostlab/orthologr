orthologr
=========

### A Software Framework for Comparative Genomics.

The `orthologr` package provides interface functions between R and common bioinformatics tools
used to perform orthology inference, multiple sequence alignments, codon alignments, dNdS estimation,
CDS prediction, and phylotranscriptomics (divergence stratigraphy).

## Fast Installation Guide

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

# install all Bioconductor packages orthologr depends on

# install Bioconductor base packages
source("http://bioconductor.org/biocLite.R")
biocLite()

# install package: Biostrings
biocLite("Biostrings")

# install package: S4Vectors
source("http://bioconductor.org/biocLite.R")
biocLite("S4Vectors")


```

## Citation


> Drost HG, Gabel A, Grosse I, Quint M. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. Mol. Biol. Evol. (2015) 32 (5): 1208-1220. doi:10.1093/molbev/msv012


## Use Cases

Learn `orthologr` by use cases provided by the following tutorials: 

- [Install Prerequisite Tools](https://github.com/HajkD/orthologr/blob/master/vignettes/Install.Rmd)
- [Perform BLAST Searches](https://github.com/HajkD/orthologr/blob/master/vignettes/blast.Rmd)
- [Perform Sequence Alignments](https://github.com/HajkD/orthologr/blob/master/vignettes/sequence_alignments.Rmd)
- [Perform Orthology Inference](https://github.com/HajkD/orthologr/blob/master/vignettes/orthology_inference.Rmd)
- [Perform dNdS Estimation](https://github.com/HajkD/orthologr/blob/master/vignettes/dNdS_estimation.Rmd)
- [Perform Divergence Stratigraphy](https://github.com/HajkD/orthologr/blob/master/vignettes/divergence_stratigraphy.Rmd)



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



