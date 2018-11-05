orthologr
=========

[![Travis-CI Build Status](https://travis-ci.org/HajkD/orthologr.svg?branch=master)](https://travis-ci.org/HajkD/orthologr)

## Comparative Genomics with R

The comparative method is a powerful approach in genomics research. Based on our knowledge about the phylogenetic relationships between species, we can study the evolution, diversification, and constraints of biological processes by comparing genomes, genes, and other genomic loci across species. The `orthologr` package aims to provide a framework to perform comparative genomics studies with R and implements R wrapper functions to the most important and most popular genomics and comparative genomics tools.

In detail, `orthologr` allows users to perform BLAST searches, orthology inference methods, multiple sequence alignments, codon alignments, dNdS estimation, and [divergence stratigraphy](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd) __between entire genomes__ with R.

The most useful implementation in `orthologr` is the ability to compute synonymous versus non-synonymous substitution rates (dN/dS)
for all orthologous genes between two __entire genomes__.


### A collection of pre-computed `dNdS maps` can be found at https://github.com/HajkD/dNdS_database

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
     delete_corrupt_cds = TRUE, # coding sequences that cannot be divided by 3 (triplets) will be removed
     ortho_detection = "RBH", # perform BLAST best reciprocal hit orthology inference
     aa_aln_type     = "pairwise", # perform pairwise global alignments of AA seqs 
     aa_aln_tool     = "NW", # using Needleman-Wunsch
     codon_aln_tool  = "pal2nal", # perform codon alignments using the tool Pal2Nal
     dnds_est.method = "Comeron", # use Comeron's method for dN/dS inference
     comp_cores      = 1 )
```

```
   query_id subject_id      dN    dS   dNdS perc_identity alig_length mismatches gap_openings
   <chr>    <chr>        <dbl> <dbl>  <dbl>         <dbl>       <dbl>      <dbl>        <dbl>
 1 AT1G010… 333554|PA… 0.106   0.254 0.420           74.0         469         80            8
 2 AT1G010… 470181|PA… 0.0402  0.104 0.388           91.1         246         22            0
 3 AT1G010… 470180|PA… 0.0150  0.126 0.118           95.5         359         12            2
 4 AT1G010… 333551|PA… 0.0135  0.116 0.116           92.0        1970         85           10
 5 AT1G010… 909874|PA… 0       0.175 0              100           213          0            0
 6 AT1G010… 470177|PA… 0.0449  0.113 0.397           89.5         648         58            5
 7 AT1G010… 918864|PA… 0.0183  0.106 0.173           95.1         366         14            2
 8 AT1G010… 909871|PA… 0.0340  0.106 0.322           90.3         300         22            2
 9 AT1G010… 470171|PA… 0.00910 0.218 0.0417          96.8         434          8            3
10 AT1G011… 333544|PA… 0.0325  0.122 0.266           93.6         528         34            0
11 AT1G011… 918858|PA… 0.00307 0.133 0.0232          99.2         529          4            0
12 AT1G011… 470161|PA… 0.00567 0.131 0.0432          98.4         453          6            1
13 AT1G011… 918855|PA… 0.13    0.203 0.641           72.6         285         68            3
14 AT1G011… 918854|PA… 0.105   0.280 0.373           84.9         179         19            2
15 AT1G011… 311317|PA… 0       0.306 0               85.6          97          0            1
16 AT1G011… 909860|PA… 0.0297  0.176 0.168           92.6         310         20            1
17 AT1G011… 311315|PA… 0.0287  0.162 0.177           94.2         533         30            1
18 AT1G012… 470156|PA… 0.0190  0.168 0.114           95.8         238         10            0
19 AT1G012… 311313|PA… 0.0207  0.154 0.134           95.3         107          5            0
20 AT1G012… 470155|PA… 0.0157  0.153 0.102           96.7        1056         35            0
# ... with 6 more variables: q_start <dbl>, q_end <dbl>, s_start <dbl>, s_end <dbl>,
#   evalue <dbl>, bit_score <dbl>
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

Users can find the corresponding map at https://github.com/HajkD/dNdS_database.

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



