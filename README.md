orthologr
=========

[![Travis-CI Build Status](https://travis-ci.org/HajkD/orthologr.svg?branch=master)](https://travis-ci.org/HajkD/orthologr)

## Comparative Genomics with R

The comparative method is a powerful approach in genomics research. Based on our knowledge about the phylogenetic relationships between species, we can study the evolution, diversification, and constraints of biological processes by comparing genomes, genes, and other genomic loci across species. The `orthologr` package aims to provide a framework to perform comparative genomics studies with R and implements R wrapper functions to the most important and most popular genomics and comparative genomics tools.

In detail, `orthologr` allows users to perform BLAST searches, orthology inference methods, multiple sequence alignments, codon alignments, dNdS estimation, and [divergence stratigraphy](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd) __between entire genomes__ with R.

The most useful implementation in `orthologr` is the ability to compute synonymous versus non-synonymous substitution rates (dN/dS)
for all orthologous genes between two __entire genomes__.


### Citation

**Please cite the following paper when using `orthologr` for your own research. This will allow me to continue
working on this software tool and will motivate me to extend its functionality and usability. Many thanks in advance :)**

> Drost et al. 2015. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ 32 (5): 1221-1231. [doi:10.1093/molbev/msv012](http://mbe.oxfordjournals.org/content/32/5/1221.abstract?sid=767aea12-1eb3-40be-8c6a-e2861f159b46)

## Install 

```r
# first install the devtools package 
install.packages("devtools")
# install orthologr on your system
source("http://bioconductor.org/biocLite.R")
biocLite("HajkD/orthologr")
```

## Use Cases

Learn `orthologr` by reading these tutorials: 

- [Install Prerequisite Tools](https://hajkd.github.io/orthologr/articles/Install.html)
- [Perform BLAST Searches](https://hajkd.github.io/orthologr/articles/blast.html)
- [Perform Sequence Alignments](https://hajkd.github.io/orthologr/articles/sequence_alignments.html)
- [Perform Orthology Inference](https://hajkd.github.io/orthologr/articles/orthology_inference.html)
- [Perform dNdS Estimation](https://hajkd.github.io/orthologr/articles/dNdS_estimation.html)
- [Perform Divergence Stratigraphy](https://hajkd.github.io/orthologr/articles/divergence_stratigraphy.html)


## Example

### Small example with internal dataset
```r
library(orthologr)

# Detect orthologous genes between a query species and a subject species
# and compute the synonymous versus non-synonymous substitution rates (dN/dS)
# following this paradigm:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) multi-core processing 'comp_cores = 1'
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1 )
```

```
       query_id            subject_id       dN     dS    dNdS
 1: AT1G01010.1 333554|PACid:16033839 0.106400 0.2537 0.41950
 2: AT1G01020.1 470181|PACid:16064328 0.040230 0.1037 0.38790
 3: AT1G01030.1 470180|PACid:16054974 0.014990 0.1265 0.11850
 4: AT1G01040.1 333551|PACid:16057793 0.013470 0.1165 0.11560
 5: AT1G01050.1 909874|PACid:16064489 0.000000 0.1750 0.00000
 6: AT1G01060.3 470177|PACid:16043374 0.044950 0.1133 0.39670
 7: AT1G01070.1 918864|PACid:16052578 0.018300 0.1059 0.17280
 8: AT1G01080.1 909871|PACid:16053217 0.033980 0.1056 0.32170
 9: AT1G01090.1 470171|PACid:16052860 0.009104 0.2181 0.04174
10: AT1G01110.2 333544|PACid:16034284 0.032480 0.1220 0.26620
11: AT1G01120.1 918858|PACid:16049140 0.003072 0.1326 0.02317
12: AT1G01140.3 470161|PACid:16036015 0.005672 0.1312 0.04324
13: AT1G01150.1 918855|PACid:16037307 0.130000 0.2028 0.64120
14: AT1G01160.1 918854|PACid:16044153 0.104600 0.2804 0.37310
15: AT1G01170.2 311317|PACid:16052302 0.000000 0.3064 0.00000
16: AT1G01180.1 909860|PACid:16056125 0.029680 0.1763 0.16830
17: AT1G01190.1 311315|PACid:16059488 0.028690 0.1618 0.17730
18: AT1G01200.1 470156|PACid:16041002 0.019050 0.1675 0.11370
19: AT1G01210.1 311313|PACid:16057125 0.020670 0.1540 0.13420
20: AT1G01220.1 470155|PACid:16047984 0.015690 0.1533 0.10230
```

### Example: Computing dN/dS values for all orthologous genes between two genomes

First, users can retrieve all coding sequences from entire genomes using the [biomartr](https://github.com/ropensci/biomartr) package ([see details here](https://ropensci.github.io/biomartr/articles/Sequence_Retrieval.html#cds-retrieval)).

```r
install.packages("biomartr")
library(biomartr)
# download all coding sequences for Mus musculus
Mmusculus_file <- biomartr::getCDS(organism = "Mus musculus", path = getwd())
# download all coding sequences for Homo sapiens
Hsapiens_file <- biomartr::getCDS(organism = "Homo sapiens", path = getwd())
# compute dN/dS values for Homo sapiens versus Mus musculus
Hs_vs_Mm_dNdS <- 
  dNdS(query_file      = Hsapiens_file,
       subject_file    = Mmusculus_file,
       ortho_detection = "RBH", 
       aa_aln_type     = "pairwise",
       aa_aln_tool     = "NW", 
       codon_aln_tool  = "pal2nal", 
       dnds_est.method = "Comeron", 
       comp_cores      = 1 )
# store result in Excel readable csv file
install.packages("readr")
readr::write_excel_csv(Hs_vs_Mm_dNdS, "Hs_vs_Mm_dNdS.csv")
```

__This way, users can compute dN/dS values for any pairwise genome comparison.__

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



