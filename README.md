orthologr
=========

## Comparative Genomics with R

### Motivation

The comparative method is a powerful approach in genomics research. Based on our knowledge about the phylogenetic relationships between species, we can study the evolution, diversification, and constraints of biological processes by comparing genomes, genes, and other genomic loci across species. The `orthologr` package aims to provide a framework to perform __large scale__ comparative genomics studies with R. `Orthologr` aims to be as easy to use as possible - from genomic data retrieval to orthology inference and dNdS estimation between several genomes.

In combination with the R package [biomartr](https://github.com/ropensci/biomartr), users can retrieve genomes, proteomes, or coding sequences for several species and use them as input
for orthology inference and dN/dS estimation with `orthologr`. The advantage of using `biomartr` in combination with `orthologr` is that
users can join the new wave of research that promotes and facilitates
[computational reproducibility in genomics studies](https://github.com/ropensci/biomartr#short-package-description) and solve the
issue of comparing genomes with different genome assembly qualities (also referred to as [genome version crisis](https://github.com/ropensci/biomartr#short-package-description)).

You can find a detailed list of all `orthologr` functions here: https://drostlab.github.io/orthologr/reference/index.html

### Citation

**Please cite the following paper in which I introduce `orthologr` when using this package for your own research. This will allow me to continue
working on this software tool and will motivate me to extend its functionality and usability in the next years. Many thanks in advance :)**

> Drost et al. 2015. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ 32 (5): 1221-1231. [doi:10.1093/molbev/msv012](http://mbe.oxfordjournals.org/content/32/5/1221.abstract?sid=767aea12-1eb3-40be-8c6a-e2861f159b46)

### Short package description

In detail, `orthologr` allows users to perform orthology inference and dN/dS estimation between two genomes or between several genomes. The following methods
to infer orthologous relationships between genes of entire genomes are available in this package:

- [__BLAST best hit__](https://drostlab.github.io/orthologr/articles/blast.html#the-blast_best-function)
- [__BLAST best reciprocal hit__](https://drostlab.github.io/orthologr/articles/blast.html#the-blast_rec-function)
- [__OrthoFinder2__](https://drostlab.github.io/orthologr/articles/orthology_inference.html)

The most useful implementation in `orthologr` is the ability to compute synonymous versus non-synonymous substitution rates (dN/dS)
for all orthologous genes between two __entire genomes__.
Available dN/dS estimation methods are:

- "NG": Nei, M. and Gojobori, T. (1986)
- "LWL": Li, W.H., et al. (1985)
- "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)
- "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)
- "YN": Yang, Z. and Nielsen, R. (2000)
- "MYN" (Modified YN): Zhang, Z., et al. (2006)
- "GMYN": Wang, D.P., et al. Biology Direct. (2009)
- "GY": Goldman, N. and Yang, Z. (1994)
- "MS": (Model Selection): based on a set of candidate models, Posada, D. (2003) 
- "MA" (Model Averaging): based on a set of candidate models, Posada, D. (2003)
- "ALL": All models toghether 

Please find more details [here](https://drostlab.github.io/orthologr/articles/dNdS_estimation.html).

## Install `orthologr`

```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# Install package dependencies
BiocManager::install(c(
        "Biostrings",
        "GenomicRanges",
        "GenomicFeatures",
        "Rsamtools",
        "rtracklayer"
))

# install metablastr from GitHub
devtools::install_github("drostlab/metablastr")

# install orthologr from GitHub
devtools::install_github("drostlab/orthologr")
```

## Use Cases

Learn `orthologr` by reading these tutorials: 

- [Install Prerequisite Tools](https://drostlab.github.io/orthologr/articles/Install.html)
- [Perform Orthology Inference and dNdS Estimation](https://drostlab.github.io/orthologr/articles/dNdS_estimation.html)
- [Perform Only Orthology Inference](https://drostlab.github.io/orthologr/articles/orthology_inference.html)
- [Perform BLAST Searches](https://drostlab.github.io/orthologr/articles/blast.html)
- [Perform Divergence Stratigraphy](https://drostlab.github.io/orthologr/articles/divergence_stratigraphy.html)
- [Perform Sequence Alignments](https://drostlab.github.io/orthologr/articles/sequence_alignments.html)



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
# A tibble: 20 x 24
   query_id subject_id      dN    dS   dNdS perc_identity num_ident_match… alig_length
   <chr>    <chr>        <dbl> <dbl>  <dbl>         <dbl>            <int>       <int>
 1 AT1G010… 333554|PA… 0.106   0.254 0.420           74.0              347         469
 2 AT1G010… 470181|PA… 0.0402  0.104 0.388           91.1              224         246
 3 AT1G010… 470180|PA… 0.0150  0.126 0.118           95.5              343         359
 4 AT1G010… 333551|PA… 0.0135  0.116 0.116           92.0             1812        1970
 5 AT1G010… 909874|PA… 0       0.175 0              100                213         213
 6 AT1G010… 470177|PA… 0.0449  0.113 0.397           89.5              580         648
 7 AT1G010… 918864|PA… 0.0183  0.106 0.173           95.1              348         366
 8 AT1G010… 909871|PA… 0.0340  0.106 0.322           90.3              271         300
 9 AT1G010… 470171|PA… 0.00910 0.218 0.0417          96.8              420         434
10 AT1G011… 333544|PA… 0.0325  0.122 0.266           93.6              494         528
11 AT1G011… 918858|PA… 0.00307 0.133 0.0232          99.2              525         529
12 AT1G011… 470161|PA… 0.00567 0.131 0.0432          98.5              446         453
13 AT1G011… 918855|PA… 0.13    0.203 0.641           72.6              207         285
14 AT1G011… 918854|PA… 0.105   0.280 0.373           84.9              152         179
15 AT1G011… 311317|PA… 0       0.306 0               85.6               83          97
16 AT1G011… 909860|PA… 0.0297  0.176 0.168           92.6              287         310
17 AT1G011… 311315|PA… 0.0287  0.162 0.177           94.2              502         533
18 AT1G012… 470156|PA… 0.0190  0.168 0.114           95.8              228         238
19 AT1G012… 311313|PA… 0.0207  0.154 0.134           95.3              102         107
20 AT1G012… 470155|PA… 0.0157  0.153 0.102           96.7             1021        1056
# … with 16 more variables: mismatches <int>, gap_openings <int>, n_gaps <int>,
#   pos_match <int>, ppos <dbl>, q_start <int>, q_end <int>, q_len <int>, qcov <int>,
#   qcovhsp <int>, s_start <int>, s_end <dbl>, s_len <dbl>, evalue <dbl>, bit_score <dbl>,
#   score_raw <dbl>
```

When running your own query file, please specify `query_file = "path/to/your/cds.fasta` instead of `system.file(..., package = "orthologr")`. The command `system.file(..., package = "orthologr")` merely references the path to the example file stored in the `orthologr` package itself.


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
       delete_corrupt_cds = FALSE,
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

Users can find the corresponding map at https://github.com/drostlab/dNdS_database.

__This way, users can compute dN/dS values for any pairwise genome comparison.__

### On Windows Systems

In some cases (when working with __WINDOWS__ machines), the installation via `devtools`
will not work properly. In this case users can try the follwing steps:

```r
# On Windows, this won't work - see ?build_github_devtools
install_github("drostlab/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first users need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards users can install devtools -> install.packages("devtools")
# and then they can run:

devtools::install_github("drostlab/orthologr", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("orthologr", lib.loc = "C:/Program Files/R/R-3.1.1/library")
```

### Troubleshooting on Windows Machines

- Install `orthologr` on a Win 8 laptop: [solution](https://github.com/drostlab/orthologr/issues/1) ( Thanks to Andres Romanowski )


## Interfaces implemented in `orthologr`:

### Perform BLAST searches with R  

* `blast()`: Perform a BLAST+ search
* `blast_best()`: Perform a BLAST+ best hit search
* `blast_rec()`: Perform a BLAST+ best reciprocal hit (BRH) search


### Perform Pairwise and Multiple Sequence Alignements with R

* `multi_aln()`: Compute Multiple Sequence Alignments based on the `clustalw`, `t_coffee`, `muscle`, `clustalo`, and `mafft` programs.
* `pairwise_aln()`: Compute Pairwise Alignments
* `codon_aln()`: Compute a Codon Alignment

### Perform Orthology Inference with R

* `orthologs()`: Main Orthology Inference Function

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

https://github.com/drostlab/orthologr/issues

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



