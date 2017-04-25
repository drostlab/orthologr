orthologr
=========

## Comparative Genomics with R

The comparative method is a powerful approach in genomics research. Based on our knowledge about the phylogenetic relationships between species, we can study the evolution, diversification, and constraints of biological processes by comparing genomes, genes, and other genomic loci across species. The `orthologr` package aims to provide a framework to perform comparative genomics studies with R and implements R wrapper functions to the most important and most popular genomics and comparative genomics tools.

In detail, `orthologr` allows users to perform BLAST searches, orthology inference methods, multiple sequence alignments, codon alignments, dNdS estimation, and [divergence stratigraphy](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd) with R.

In future versions we will include a vast range of analytics functions to perform state-of-the art and cutting edge comparative genomics studies.

### Citation

**Please cite the following paper when using `orthologr` for your own research. This will allow me to continue
working on this software tool and will motivate me to extend its functionality and usability. Many thanks in advance :)**

> Drost HG, Gabel A, Grosse I, Quint M. 2015. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ 32 (5): 1221-1231. [doi:10.1093/molbev/msv012](http://mbe.oxfordjournals.org/content/32/5/1221.abstract?sid=767aea12-1eb3-40be-8c6a-e2861f159b46)

## Use Cases

Learn `orthologr` by reading these tutorials: 

- [Install Prerequisite Tools](https://github.com/HajkD/orthologr/blob/master/vignettes/Install.Rmd)
- [Perform BLAST Searches](https://github.com/HajkD/orthologr/blob/master/vignettes/blast.Rmd)
- [Perform Sequence Alignments](https://github.com/HajkD/orthologr/blob/master/vignettes/sequence_alignments.Rmd)
- [Perform Orthology Inference](https://github.com/HajkD/orthologr/blob/master/vignettes/orthology_inference.Rmd)
- [Perform dNdS Estimation](https://github.com/HajkD/orthologr/blob/master/vignettes/dNdS_estimation.Rmd)
- [Perform Divergence Stratigraphy](https://github.com/HajkD/orthologr/blob/master/vignettes/divergence_stratigraphy.Rmd)


## Install

```r
# first install the devtools package 
install.packages("devtools")
# install orthologr on your system
source("http://bioconductor.org/biocLite.R")
biocLite("HajkD/orthologr")
```

### On Windows Systems

In some cases (when working with __WINDOWS__ machines), the installation via `devtools`
will not work properly. In this case users can try the follwing steps:

```r
# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first users need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards users can install devtools -> install.packages("devtools")
# and then they can run:

devtools::install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("orthologr", lib.loc = "C:/Program Files/R/R-3.1.1/library")
```

### Troubleshooting on Windows Machines

- Install `orthologr` on a Win 8 laptop: [solution](https://github.com/HajkD/orthologr/issues/1) ( Thanks to Andres Romanowski )


## Interfaces implemented in `orthologr`:

### Perform BLAST searches with R  

* `advanced_blast()`: Perform an advanced BLAST+ search
* `advanced_makedb()`: Create a BLASTable database with makeblastdb (advanced options)
* `blast()`: Perform a BLAST+ search
* `blast.nr()`: Perform a BLASTP search against NCBI nr
* `blast_best()`: Perform a BLAST+ best hit search
* `blast_rec()`: Perform a BLAST+ best reciprocal hit (BRH) search
* `delta.blast()`: Perform a DELTA-BLAST Search


### Perform Pairwise and Multiple Sequence Alignements with R

* `multi_aln()`: Compute Multiple Sequence Alignments based on the `clustalw`, `t_coffee`, `muscle`, `clustalo`, and `mafft` programs.
* `pairwise_aln()`: Compute Pairwise Alignments
* `codon_aln()`: Compute a Codon Alignment

### Perform Orthology Inference with R

* `orthologs()`: Main Orthology Inference Function
* `ProteinOrtho()`: Orthology Inference with ProteinOrtho

### Perform Population Genomics with R

* `dNdS()`: Compute dNdS values for two organisms
* `substitutionrate()`: Internal function for dNdS computations

### Read and Write CDS, Genomes, and Proteomes

* `read.cds()`: Read the CDS of a given organism
* `read.genome()`: Read the genome of a given organism
* `read.proteome()`: Read the proteome of a given organism
* `write.proteome()`: Save a proteome in fasta format


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/HajkD/orthologr/issues

## Licenses

The `orthologr` package includes source code that has been published under following licenses:

### gestimator

```
All files included in `orthologr` that were taken from gestimator are 
also part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
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
version 2 of the License, or any later version.

```

### KaKs_Calculator

In `orthologr` the file `parseFastaIntoAXT.pl` is stored in `/inst/KaKs_Calc_parser`.

```
The parseFastaIntoAXT.pl script is freely available under GNU GPL v3 
Licence and included in the KaKs_Calculator project that can be found at 
https://code.google.com/p/kaks-calculator/

```



