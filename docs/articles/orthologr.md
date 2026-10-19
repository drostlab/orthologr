# orthologr

`orthologr` infers orthologous genes between genomes, estimates their
dN/dS, and clusters proteins into homologous groups, all without leaving
the R environment. Four steps get you started.

## 1. Install `orthologr`

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "GenomicFeatures",
                       "Rsamtools", "rtracklayer", "pwalign"))
install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools"))
devtools::install_github("drostlab/metablastr")
devtools::install_github("drostlab/orthologr")
```

Windows users also need
[Rtools](https://drostlab.github.io/orthologr/articles/setup_windows.html#install-rtools).

## 2. Install DIAMOND2

``` r

library(orthologr)
diamond_install()
```

This downloads `DIAMOND2` for your operating system into a user cache
folder, with no administrator rights or `PATH` changes. It is the only
external tool needed for steps 3 and 4. See
[Install](https://drostlab.github.io/orthologr/articles/Install.html#diamond-install)
for alternatives (conda, manual install, BLAST+).

## 3. Estimate dN/dS with `dNdS()`

[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) finds
orthologous gene pairs between two species and computes dN, dS and dN/dS
for each pair from their coding sequences.

``` r

dnds <- dNdS(
  query_file     = system.file("seqs/ortho_thal_cds.fasta", package = "orthologr"),
  subject_file   = system.file("seqs/ortho_lyra_cds.fasta", package = "orthologr"),
  print_citation = FALSE
)
#> 
#> Starting orthology inference (RBH) and dNdS estimation (Comeron) using the follwing parameters:
#> query = 'ortho_thal_cds.fasta'
#> subject = 'ortho_lyra_cds.fasta'
#> aligner = 'diamond'
#> sensitivity_mode = 'fast'
#> seq_type = 'cds'
#> e-value: 1E-5
#> aa_aln_type = 'pairwise'
#> aa_aln_tool = 'NW'
#> comp_cores = '1'
#> 
#> Starting Orthology Inference ...
#> Running diamond version 2.2.1 ...
#> sensitivity mode: fast
#> creating a diamond database
#> Starting DIAMOND2 search ...
#> DIAMOND2 search completed in 0.03 sec.
#> Running diamond version 2.2.1 ...
#> sensitivity mode: fast
#> creating a diamond database
#> Starting DIAMOND2 search ...
#> DIAMOND2 search completed in 0.03 sec.
#> Orthology Inference Completed.
#> Starting dN/dS Estimation ...
#> dN/dS Estimation Completed.

head(dnds[, c("query_id", "subject_id", "dN", "dS", "dNdS")])
#> # A tibble: 6 × 5
#>   query_id    subject_id                dN    dS  dNdS
#>   <chr>       <chr>                  <dbl> <dbl> <dbl>
#> 1 AT1G01010.1 333554|PACid:16033839 0.106  0.254 0.420
#> 2 AT1G01020.1 470181|PACid:16064328 0.0402 0.104 0.388
#> 3 AT1G01030.1 470180|PACid:16054974 0.0150 0.126 0.118
#> 4 AT1G01040.1 333551|PACid:16057793 0.0135 0.116 0.116
#> 5 AT1G01050.1 909874|PACid:16064489 0      0.175 0    
#> 6 AT1G01060.3 470177|PACid:16043374 0.0449 0.113 0.397
```

More: [Orthology
inference](https://drostlab.github.io/orthologr/articles/orthology_inference.md),
[dNdS
estimation](https://drostlab.github.io/orthologr/articles/dNdS_estimation.md),
[Divergence
stratigraphy](https://drostlab.github.io/orthologr/articles/divergence_stratigraphy.md).

## 4. Cluster proteins with `deepclust()`

[`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md)
groups the proteins of one or more proteomes into homologous clusters.

``` r

clusters <- deepclust(
  input_file = c(
    system.file("seqs/ortho_thal_aa.fasta", package = "orthologr"),
    system.file("seqs/ortho_lyra_aa.fasta", package = "orthologr")
  ),
  mutual_cover = 70 # 70% mutual coverage of the aligned sequences
)
#> Running diamond version 2.2.1 ...
#> Merging 2 input files ...
#> Merged 40 sequences into: deepclust_merged_2_files.fasta
#> Building DIAMOND2 database for deepclust ...
#> Running diamond deepclust ...
#> diamond deepclust completed in 0.77 sec.

head(clusters)
#> # A tibble: 6 × 2
#>   representative_id     member_id            
#>   <chr>                 <chr>                
#> 1 333551|PACid:16057793 333551|PACid:16057793
#> 2 333551|PACid:16057793 AT1G01040.1          
#> 3 470155|PACid:16047984 AT1G01220.1          
#> 4 470155|PACid:16047984 470155|PACid:16047984
#> 5 470177|PACid:16043374 AT1G01060.3          
#> 6 470177|PACid:16043374 470177|PACid:16043374
```

More: [Protein clustering with
DeepClust](https://drostlab.github.io/orthologr/articles/deepclust.md).
