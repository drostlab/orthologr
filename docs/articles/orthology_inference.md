# Orthology Inference using orthologr

Orthology Inference is the process of detecting [orthologous
genes](https://www.biostars.org/p/1595/) between a query organism and a
set of subject organisms. Motivated by the discussion on which orthology
inference method, paradigm or program is the [most
accurate](https://www.biostars.org/p/7568/) or [the
fastest](https://www.nature.com/articles/nmeth.3830) one, the
`orthologr` package provides functions to perform genome wide orthology
inference following two different paradigms:

> 1.  orthology inference based on 1-1 orthology relationship using
>     `DIAMOND2 reciprocal best hit` (default, recommended) or
>     `BLAST reciprocal best hit`
>
> 2.  orthology inference based on 1-many homology relationship using
>     orthologous groups (OrthoFinder2)

To perform orthology inference you can start with the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function provided by `orthologr`. The
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function takes nucleotide or protein sequences stored in fasta files for
a set of organisms as input and performs orthology inference to detect
orthologous genes or orthologous groups within the given organisms based
on selected orthology inference programs.

The following interfaces are implemented in the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function:

### DIAMOND2 based methods (recommended):

- DIAMOND2 reciprocal best hit (DIAMOND_RBH): see
  [`diamond_rec()`](https://drostlab.github.io/orthologr/reference/diamond_rec.md).
  **This is the default** (`ortho_detection = "DIAMOND_RBH"`). DIAMOND2
  is up to 10,000× faster than BLAST in default mode and retains full
  BLAST-level sensitivity in `ultra-sensitive` mode.

- DIAMOND2 best hit (DIAMOND_BH): see
  [`diamond_best()`](https://drostlab.github.io/orthologr/reference/diamond_best.md).

### BLAST based methods:

- BLASTp best hit (BH): see
  [`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md).

- BLASTp reciprocal best hit (RBH): see
  [`blast_rec()`](https://drostlab.github.io/orthologr/reference/blast_rec.md).

### Orthogroup based methods:

- [Orthofinder2](https://github.com/davidemms/OrthoFinder): see
  [`orthofinder2()`](https://drostlab.github.io/orthologr/reference/orthofinder2.md).

### Clustering based methods:

- [DIAMOND DeepClust](https://github.com/bbuchfink/diamond): ultra-fast
  protein sequence clustering that retains sensitivity at low sequence
  identity. Useful for large proteome comparisons where pairwise hit
  methods become expensive. See
  [`deepclust()`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  and the [DeepClust
  vignette](https://drostlab.github.io/orthologr/articles/deepclust.html)
  for details.

### Examples: DIAMOND2 based methods (default)

Using a simple example stored in the package environment of `orthologr`
you can get an impression on how to use the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function.

**Note:** it is assumed that when using
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
all corresponding programs you want to use are already installed on your
machine and are executable via either the default execution PATH or you
specifically define the location of the executable file via the `path`
argument that can be passed to
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md).
For DIAMOND2-based methods (the default), make sure
[DIAMOND2](https://github.com/bbuchfink/diamond) is installed. See the
[Installation
Vignette](https://drostlab.github.io/orthologr/articles/Install.html)
for details.

In the following examples, we will use toy fasta files stored within the
`orthologr` package:
`system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr')`.

``` r

# perform orthology inference using DIAMOND2 reciprocal best hit (default)
# and fasta sequence files storing protein sequences
orthologr::orthologs(query_file    = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          comp_cores      = 1,
          clean_folders   = FALSE)
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
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974          93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793          93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374          87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578          92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217          89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284          87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153          78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984          96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

This small example returns 20 orthologous genes between *Arabidopsis
thaliana* and *Arabidopsis lyrata*. As you can see, the `query_file` and
`subject_files` arguments take the proteomes of *Arabidopsis thaliana*
(`query_file`) and *Arabidopsis lyrata* (`subject_files`) stored in
fasta files. The `seq_type` argument specifies that you will pass
protein sequences (proteomes) to the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function. In case you only have either genomes (DNA sequences) or CDS
files, you can modify the `seq_type` argument to `seq_type = "dna"`
(when working with only genome data) or `seq_type = "cds"` (when working
with CDS files). Internally the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function will perform a CDS prediction using `predict_cds()` and will
furthermore translate the predicted CDS sequences into protein
sequences. Analogously when `seq_type = "cds"` is specified, internally
the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function will translate all CDS sequences into protein sequences to run
orthology inference based on protein sequences.

**The advantage of this type of output is that in addition to having the
orthology relationship between genes from two different genomes, users
also retrieve the DIAMOND2 alignment results of the respective orthology
relationship (in the same format as BLAST) which allows them to perform
subsequent filtering for either more conservative or more liberal
orthology relationships.**

**Note**: future versions of `orthologr` will allow to perform orthology
inference using DNA sequences. Nevertheless, since most orthology
inference methods or paradigms rely on protein sequences, the first
version of `orthologr` will follow this paradigm.

In case you have to specify the path to DIAMOND2 you can use the `path`
argument as follows:

``` r

# using an external execution path for DIAMOND2
orthologr::orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          path            = "here/path/to/diamond",
          clean_folders   = FALSE,
          comp_cores      = 1)
```

When you are working on a multicore machine, you can also specify the
`comp_cores` argument that will allow you to run all analyses in
parallel (to speed up computations).

``` r

# running DIAMOND2 orthology inference in parallel using 2 cores
orthologr::orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          clean_folders   = FALSE,
          comp_cores      = 2)
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
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974          93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793          93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374          87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578          92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217          89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284          87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153          78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984          96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

In this case 2 cores are being used to perform parallel processing, the
`clean_folders` argument specifies that all files returned by the
corresponding orthology inference method are removed after analyses.

## Program specific use of the orthologs() function

In this section small examples will illustrate the use of the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function for each orthology inference program.

### DIAMOND2 best hit

The DIAMOND2 best hit method is a uni-directional DIAMOND2 best hit
search of a `query organism A` against a `subject organism B` based on
the `e-value`. It uses the same hit-filtering logic as
[`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md)
but is substantially faster.

Orthology Inference using DIAMOND2 best hit can be performed by
specifying the argument `ortho_detection = "DIAMOND_BH"` and one
computing core `comp_core = 1`:

``` r

# perform orthology inference using DIAMOND2 best hit
# and fasta sequence files storing protein sequences
orthologr::orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_BH",
          clean_folders   = TRUE,
          comp_cores      = 1)
#> Running diamond version 2.2.1 ...
#> sensitivity mode: fast
#> creating a diamond database
#> Starting DIAMOND2 search ...
#> DIAMOND2 search completed in 0.03 sec.
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974          93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793          93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374          87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578          92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217          89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284          87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153          78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984          96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

The resulting table stores the orthologous gene pairs and the
corresponding alignment statistics. DIAMOND2 returns additional columns
compared to BLAST (e.g. `n_gaps`, `pos_match`, `ppos`, `q_len`, `qcov`,
`s_len`, `score_raw`).

### DIAMOND2 reciprocal best hit (default, recommended)

The DIAMOND2 reciprocal best hit is a bi-directional DIAMOND2 best hit
search of a `query organism A` against a `subject organism B` based on
the `e-value`. This is the **default** method in
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md).

The algorithm runs as follows:

1.  `bh_A <- diamond_best(A, B)`

2.  `bh_B <- diamond_best(B, A)`

3.  `join(bh_A, bh_B)` by tuple `(query_id, subject_id)`

Only when the tuple `(query_id, subject_id)` is returned as best hit in
both DIAMOND2 directions is it retained as an orthologous gene pair.

``` r

library(orthologr)

# perform orthology inference using DIAMOND2 reciprocal best hit (default)
# and fasta sequence files storing protein sequences
orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          clean_folders   = FALSE,
          comp_cores      = 1)
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
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974          93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793          93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374          87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578          92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217          89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284          87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153          78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984          96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

You can also adjust the DIAMOND2 sensitivity level using the
`sensitivity_mode` argument. The default is `"fast"`, but for more
distant species comparisons you may want `"sensitive"` or
`"ultra-sensitive"`:

``` r

library(orthologr)

# DIAMOND2 RBH with ultra-sensitive mode (equivalent to BLAST sensitivity)
orthologs(query_file       = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files    = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type         = "protein",
          ortho_detection  = "DIAMOND_RBH",
          sensitivity_mode = "ultra-sensitive",
          clean_folders    = FALSE,
          comp_cores       = 1)
#> Running diamond version 2.2.1 ...
#> sensitivity mode: ultra-sensitive
#> creating a diamond database
#> Starting DIAMOND2 search ...
#> DIAMOND2 search completed in 0.33 sec.
#> Running diamond version 2.2.1 ...
#> sensitivity mode: ultra-sensitive
#> creating a diamond database
#> Starting DIAMOND2 search ...
#> DIAMOND2 search completed in 0.33 sec.
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974          93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793          93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374          87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578          92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217          89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284          87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153          78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984          96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

The full set of sensitivity modes available in DIAMOND2 is:

| `sensitivity_mode` | Speed | Use case |
|----|----|----|
| `"faster"` | Fastest | Hits \>70% identity |
| `"fast"` | Fast (default) | Hits \>70% identity |
| `"mid-sensitive"` | Moderate | Between fast and sensitive |
| `"sensitive"` | Sensitive | Full sensitivity for hits \>40% identity |
| `"more-sensitive"` | More sensitive | — |
| `"very-sensitive"` | Very sensitive | — |
| `"ultra-sensitive"` | BLAST-equivalent | Most sensitive; slowest |

To retrieve the full DIAMOND2 hit table for subsequent analyses, set
`clean_folders = FALSE`. The corresponding hit table can then be found
in `file.path(tempdir(), "_blast_db")`.

``` r

library(orthologr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# perform orthology inference using DIAMOND2 reciprocal best hit
RBH <- orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          clean_folders   = FALSE,
          comp_cores      = 1)
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

glimpse(RBH)
#> Rows: 20
#> Columns: 21
#> Groups: query_id [20]
#> $ query_id          <chr> "AT1G01010.1", "AT1G01020.1", "AT1G01030.1", "AT1G01040.1", "AT1G01050.1", "AT1G01060.3", "AT1G01070.1", "AT1G01080.1", "AT1G01090.1", "AT1G01110.2", "AT1G01120.1", "AT1G01140.3", "AT1G01150.1", "AT1G01160.1", "AT1G01170.2", "AT1G01180.1", "AT1G01190.1", "AT1G01200.1", "AT1G01210.1", "AT1G01220.1"
#> $ subject_id        <chr> "333554|PACid:16033839", "470181|PACid:16064328", "470180|PACid:16054974", "333551|PACid:16057793", "909874|PACid:16064489", "470177|PACid:16043374", "918864|PACid:16052578", "909871|PACid:16053217", "470171|PACid:16052860", "333544|PACid:16034284", "918858|PACid:16049140", "470161|PACid:16036015", "918855|PACid:16037307", "918854|PACid:16044153", "311317|PACid:16052302", "909860|PACid:16056125", "311315|PACid:16059488", "470156|PACid:16041002", "311313|PACid:16057125", "470155|PACid:16047984"
#> $ perc_identity     <dbl> 73.2, 91.1, 93.3, 93.4, 100.0, 87.5, 92.6, 89.3, 96.8, 87.7, 99.2, 98.5, 72.6, 78.8, 85.6, 92.6, 94.2, 95.8, 95.3, 96.5
#> $ num_ident_matches <int> 347, 224, 335, 1840, 213, 567, 339, 268, 420, 463, 525, 446, 207, 141, 83, 287, 502, 228, 102, 1019
#> $ alig_length       <int> 474, 246, 359, 1969, 213, 648, 366, 300, 434, 528, 529, 453, 285, 179, 97, 310, 533, 238, 107, 1056
#> $ mismatches        <int> 75, 22, 20, 58, 0, 71, 23, 25, 8, 65, 4, 6, 68, 30, 0, 20, 30, 10, 5, 37
#> $ gap_openings      <int> 8, 0, 2, 7, 0, 5, 2, 2, 3, 0, 0, 1, 3, 2, 1, 1, 1, 0, 0, 0
#> $ n_gaps            <int> 52, 0, 4, 71, 0, 10, 4, 7, 6, 0, 0, 1, 10, 8, 14, 3, 1, 0, 0, 0
#> $ pos_match         <int> 369, 231, 338, 1870, 213, 586, 342, 275, 425, 475, 527, 448, 234, 145, 83, 296, 518, 233, 105, 1035
#> $ ppos              <dbl> 77.8, 93.9, 94.2, 95.0, 100.0, 90.4, 93.4, 91.7, 97.9, 90.0, 99.6, 98.9, 82.1, 81.0, 85.6, 95.5, 97.2, 97.9, 98.1, 98.0
#> $ q_start           <int> 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 5, 4, 2, 1, 1, 1, 1, 1
#> $ q_end             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 429, 528, 529, 452, 289, 175, 84, 307, 533, 238, 107, 1056
#> $ q_len             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 429, 528, 529, 452, 346, 196, 84, 479, 536, 238, 107, 1056
#> $ qcov              <dbl> 100.0, 100.0, 100.0, 99.7, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 82.4, 87.8, 98.8, 64.1, 99.4, 100.0, 100.0, 100.0
#> $ qcovhsp           <dbl> 100.0, 100.0, 100.0, 99.7, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 82.4, 87.8, 98.8, 64.1, 99.4, 100.0, 100.0, 100.0
#> $ s_start           <int> 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 16, 2, 4, 1, 6, 1, 1, 1
#> $ s_end             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 433, 528, 529, 453, 290, 179, 100, 310, 537, 238, 107, 1056
#> $ s_len             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 433, 528, 529, 453, 294, 225, 100, 316, 540, 238, 107, 1056
#> $ evalue            <dbl> 1.37e-211, 2.27e-156, 4.01e-210, 0.00e+00, 4.47e-158, 0.00e+00, 6.01e-227, 1.14e-181, 5.41e-305, 1.18e-309, 0.00e+00, 0.00e+00, 5.31e-148, 1.62e-94, 8.03e-57, 4.57e-214, 0.00e+00, 6.96e-163, 2.04e-80, 0.00e+00
#> $ bit_score         <dbl> 582, 426, 571, 3544, 427, 1028, 614, 494, 817, 837, 1029, 865, 410, 266, 161, 584, 1007, 441, 222, 2016
#> $ score_raw         <dbl> 1500, 1094, 1471, 9190, 1098, 2658, 1583, 1271, 2110, 2162, 2661, 2235, 1055, 680, 408, 1506, 2603, 1135, 566, 5224
```

### BLASTp best hit

The BLASTp best hit method is a uni-directional BLAST best hit search of
a `query organism A` against a `subject organism B` based on the
`e-value`.

Orthology Inference using BLASTp best hit can be performed by specifying
the argument `ortho_detection = "BH"` and one computing core
`comp_core = 1`:

``` r

# perform orthology inference using BLAST best hit
# and fasta sequence files storing protein sequences
orthologr::orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein", 
          ortho_detection = "BH", 
          clean_folders   = TRUE,
          comp_cores      = 1)
#> Running blastp: 2.17.0+ ...
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          74.0               347         469         80            8     42       370  78.9       1   430   430   100     100       1   466   466 0               627      1617
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246   100     100       1   246   246 8.91e-168       454      1169
#>  3 AT1G01030.1 470180|PACid:16054974          95.5               343         359         12            2      4       346  96.4       1   359   359   100     100       1   355   355 0               698      1801
#>  4 AT1G01040.1 333551|PACid:16057793          92.0              1812        1970         85           10     73      1853  94.1       6  1910  1910    99      99       2  1963  1963 0              3704      9604
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213   100     100       1   213   213 3.27e-162       437      1125
#>  6 AT1G01060.3 470177|PACid:16043374          89.5               580         648         58            5     10       601  92.8       1   646   646   100     100       1   640   640 0              1037      2681
#>  7 AT1G01070.1 918864|PACid:16052578          95.1               348         366         14            2      4       351  95.9       1   366   366   100     100       1   362   362 0               696      1797
#>  8 AT1G01080.1 909871|PACid:16053217          90.3               271         300         22            2      7       279  93         1   294   294   100     100       1   299   299 1.21e-180       491      1264
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429   100     100       1   433   433 0               859      2220
#> 10 AT1G01110.2 333544|PACid:16034284          93.6               494         528         34            0      0       507  96.0       1   528   528   100     100       1   528   528 0               972      2512
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529   100     100       1   529   529 0              1092      2825
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452   100     100       1   453   453 0               918      2372
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346    82      82      16   290   294 4.13e-152       421      1082
#> 14 AT1G01160.1 918854|PACid:16044153          84.9               152         179         19            2      8       157  87.7       4   175   196    88      88       2   179   225 2.39e- 95       268       685
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84    99      99       4   100   100 1.18e- 55       158       400
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479    64      64       1   310   316 0               576      1484
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536    99      99       6   537   540 0              1036      2679
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238   100     100       1   238   238 3.85e-174       470      1209
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107   100     100       1   107   107 1.81e- 77       215       547
#> 20 AT1G01220.1 470155|PACid:16047984          96.7              1021        1056         35            0      0      1037  98.2       1  1056  1056   100     100       1  1056  1056 0              2106      5456
```

The resulting table stores the orthologous gene pairs and the
corresponding BLAST alignment statistics.

### BLASTp best reciprocal hit

The BLAST best reciprocal hit is a bi-directional BLAST best hit search
of a `query organism A` against a `subject organism B` based on the
`e-value`.

The Algorithm for BLAST best reciprocal hit runs as follows:

1.  `bh_A <- best_hit(A,B)`

2.  `bh_B <- best_hit(B,A)`

3.  `join(bh_A,bh_B)` by tupel `(query_id, subject_id)`

In other words, only in case the tuple `(query_id, subject_id)` is
returned as best hit based on the `e-value` in both BLAST directions,
the corresponding tupel `(query_id, subject_id)` is retained as
orthologous gene pair.

``` r


library(orthologr)

# perform orthology inference using BLAST reciprocal best hit
# and fasta sequence files storing protein sequences
orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein", 
          ortho_detection = "RBH", 
          clean_folders   = FALSE,
          comp_cores      = 1)
#> Running blastp: 2.17.0+ ...
#> Running blastp: 2.17.0+ ...
#> # A tibble: 20 × 21
#> # Groups:   query_id [20]
#>    query_id    subject_id            perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839          74.0               347         469         80            8     42       370  78.9       1   430   430   100     100       1   466   466 0               627      1617
#>  2 AT1G01020.1 470181|PACid:16064328          91.1               224         246         22            0      0       231  93.9       1   246   246   100     100       1   246   246 8.91e-168       454      1169
#>  3 AT1G01030.1 470180|PACid:16054974          95.5               343         359         12            2      4       346  96.4       1   359   359   100     100       1   355   355 0               698      1801
#>  4 AT1G01040.1 333551|PACid:16057793          92.0              1812        1970         85           10     73      1853  94.1       6  1910  1910    99      99       2  1963  1963 0              3704      9604
#>  5 AT1G01050.1 909874|PACid:16064489         100                 213         213          0            0      0       213 100         1   213   213   100     100       1   213   213 3.27e-162       437      1125
#>  6 AT1G01060.3 470177|PACid:16043374          89.5               580         648         58            5     10       601  92.8       1   646   646   100     100       1   640   640 0              1037      2681
#>  7 AT1G01070.1 918864|PACid:16052578          95.1               348         366         14            2      4       351  95.9       1   366   366   100     100       1   362   362 0               696      1797
#>  8 AT1G01080.1 909871|PACid:16053217          90.3               271         300         22            2      7       279  93         1   294   294   100     100       1   299   299 1.21e-180       491      1264
#>  9 AT1G01090.1 470171|PACid:16052860          96.8               420         434          8            3      6       425  97.9       1   429   429   100     100       1   433   433 0               859      2220
#> 10 AT1G01110.2 333544|PACid:16034284          93.6               494         528         34            0      0       507  96.0       1   528   528   100     100       1   528   528 0               972      2512
#> 11 AT1G01120.1 918858|PACid:16049140          99.2               525         529          4            0      0       527  99.6       1   529   529   100     100       1   529   529 0              1092      2825
#> 12 AT1G01140.3 470161|PACid:16036015          98.5               446         453          6            1      1       448  98.9       1   452   452   100     100       1   453   453 0               918      2372
#> 13 AT1G01150.1 918855|PACid:16037307          72.6               207         285         68            3     10       234  82.1       5   289   346    82      82      16   290   294 4.13e-152       421      1082
#> 14 AT1G01160.1 918854|PACid:16044153          84.9               152         179         19            2      8       157  87.7       4   175   196    88      88       2   179   225 2.39e- 95       268       685
#> 15 AT1G01170.2 311317|PACid:16052302          85.6                83          97          0            1     14        83  85.6       2    84    84    99      99       4   100   100 1.18e- 55       158       400
#> 16 AT1G01180.1 909860|PACid:16056125          92.6               287         310         20            1      3       296  95.5       1   307   479    64      64       1   310   316 0               576      1484
#> 17 AT1G01190.1 311315|PACid:16059488          94.2               502         533         30            1      1       518  97.2       1   533   536    99      99       6   537   540 0              1036      2679
#> 18 AT1G01200.1 470156|PACid:16041002          95.8               228         238         10            0      0       233  97.9       1   238   238   100     100       1   238   238 3.85e-174       470      1209
#> 19 AT1G01210.1 311313|PACid:16057125          95.3               102         107          5            0      0       105  98.1       1   107   107   100     100       1   107   107 1.81e- 77       215       547
#> 20 AT1G01220.1 470155|PACid:16047984          96.7              1021        1056         35            0      0      1037  98.2       1  1056  1056   100     100       1  1056  1056 0              2106      5456
```

In case you would like to store the corresponding `hit tables` returned
by BLAST for subsequent analyses, you can specify the
`clean_folders = FALSE` argument. The corresponding BLAST hit table can
then be found in `file.path(tempdir(),"_blast_db")`.

A detailed overview of further analyses that can be done with the
corresponding BLAST output can be found in the [BLAST
vignette](https://drostlab.github.io/orthologr/articles/blast.html).
