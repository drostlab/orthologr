orthologr
=========

### A package combining common gene orthology detection methods and dNdS models for pipeline processing 

## Installing prerequisite tools

The `orthologr` package provides interfaces between R and common bioinformatics tools
used to perform orthology inference, multiple sequence alignments, codon alignments, dNdS computations,
and CDS predictions.


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