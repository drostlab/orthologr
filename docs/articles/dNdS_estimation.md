# Orthology inference dNdS estimation of orthologous genes using orthologr

The dN/dS ratio quantifies the mode and strength of selection acting on
a pair of orthologous genes. This selection pressure can be quantified
by comparing synonymous substitution rates (dS) that are assumed to be
neutral with nonsynonymous substitution rates (dN), which are exposed to
selection as they change the amino acid composition of a protein ([Mugal
et al., 2013](https://academic.oup.com/mbe/article/31/1/212/1049986)).

The `orthologr` package provides a function named
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) to
perform dNdS estimation on pairs of orthologous genes. The
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
function takes the CDS files of two organisms of interest (`query_file`
and `subject_file`) and computes the dNdS estimation values for
orthologous gene pairs between these organisms.

**Note:** the following dNdS estimation methods are based on
KaKs_Calculator 2:

- “NG”: Nei, M. and Gojobori, T. (1986)

- “LWL”: Li, W.H., et al. (1985)

- “LPB”: Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)

- “MLWL” (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)

- “YN”: Yang, Z. and Nielsen, R. (2000)

- “MYN” (Modified YN): Zhang, Z., et al. (2006)

- “GMYN”: Wang, D.P., et al. Biology Direct. (2009)

- “GY”: Goldman, N. and Yang, Z. (1994)

- “MS”: (Model Selection): based on a set of candidate models,
  Posada, D. (2003)

- “MA” (Model Averaging): based on a set of candidate models, Posada, D.
  (2003)

- “ALL”: All models toghether

It is assumed that when you choose one of these dNdS estimation methods
you have [KaKs_Calculator
2](https://sourceforge.net/projects/kakscalculator2/) installed on your
machine and it can be executed from the default execution `PATH`.

The following pipeline resembles an example dNdS estimation procedure:

1.  Orthology Inference: e.g. DIAMOND2 reciprocal best hit (RBH) —
    default; or BLAST reciprocal best hit

2.  Pairwise sequence alignment: e.g. clustalw for pairwise amino acid
    sequence alignments

3.  Codon Alignment: e.g. pal2nal program

4.  dNdS estimation: e.g. [Yang, Z. and Nielsen, R.
    (2000)](https://academic.oup.com/mbe/article/17/1/32/975527) (YN)

**Note:** it is assumed that when using
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) all
corresponding programs you want to use are already installed on your
machine and are executable via either the default execution PATH or you
specifically define the location of the executable program via the
`aa_aln_path` or `aligner_path` argument that can be passed to
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md). By
default
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) uses
DIAMOND2 (`aligner = "diamond"`); you can switch to BLAST by setting
`aligner = "blast"`. See the [Sequence Alignments
vignette](https://drostlab.github.io/orthologr/articles/sequence_alignments.html)
for details.

The following example shall illustrate a dNdS estimation process.

``` r

library(orthologr)
# get a dNdS table using:
# 1) DIAMOND2 reciprocal best hit for orthology inference (default; aligner = "diamond", ortho_detection = "RBH")
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single core processing 'comp_cores = 1'
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     aligner         = "diamond",
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1)
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
#> 
#> Please cite the following paper when using orthologr for your own research:
#> Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. 2015. Mol. Biol. Evol. 32 (5): 1221-1231.
#> 
#> 
#> # A tibble: 20 × 24
#>    query_id    subject_id                 dN    dS   dNdS perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                   <dbl> <dbl>  <dbl>         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839 0.106   0.254 0.420           73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328 0.0402  0.104 0.388           91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974 0.0150  0.126 0.118           93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793 0.0135  0.116 0.116           93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489 0       0.175 0              100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374 0.0449  0.113 0.397           87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578 0.0183  0.106 0.173           92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217 0.0340  0.106 0.322           89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860 0.00910 0.218 0.0417          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284 0.0325  0.122 0.266           87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140 0.00307 0.133 0.0232          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015 0.00567 0.131 0.0432          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307 0.13    0.203 0.641           72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153 0.105   0.280 0.373           78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302 0       0.306 0               85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125 0.0297  0.176 0.168           92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488 0.0287  0.162 0.177           94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002 0.0190  0.168 0.114           95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125 0.0207  0.154 0.134           95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984 0.0157  0.153 0.102           96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

Some outputs include `NA` values. To filter for `NA` values or a
specific `dnds.threshold`, you can use the
[`filter_dNdS()`](https://drostlab.github.io/orthologr/reference/filter_dNdS.md)
function. The
[`filter_dNdS()`](https://drostlab.github.io/orthologr/reference/filter_dNdS.md)
function takes the output data.table returned by
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) and
filters the output by the following criteria:

1.  all dN values having an NA value are omitted

2.  all dS values having an NA value are omitted

3.  all dNdS values \>= the specified `dnds.threshold` are omitted

``` r

library(orthologr)
# get dNdS estimated for orthologous genes between A. thaliana and A. lyrata
Ath_Aly_dnds <- dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     aligner         = "diamond",
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1)
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
#> 
#> Please cite the following paper when using orthologr for your own research:
#> Drost et al. Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis. 2015. Mol. Biol. Evol. 32 (5): 1221-1231.
#> 
#> 
# filter for:
# 1) all dN values having an NA value are omitted
# 2) all dS values having an NA value are omitted
# 3) all dNdS values >= 2 are omitted
filter_dNdS(Ath_Aly_dnds, dnds.threshold = 2)
#> Filtering out NA values in dN or dS and all values with dNdS > 2 ...
#> Initial input contains 20 rows.
#> Filtering done. New output table contains 20 rows.
#> # A tibble: 20 × 24
#>    query_id    subject_id                 dN    dS   dNdS perc_identity num_ident_matches alig_length mismatches gap_openings n_gaps pos_match  ppos q_start q_end q_len  qcov qcovhsp s_start s_end s_len    evalue bit_score score_raw
#>    <chr>       <chr>                   <dbl> <dbl>  <dbl>         <dbl>             <int>       <int>      <int>        <int>  <int>     <int> <dbl>   <int> <int> <int> <dbl>   <dbl>   <int> <int> <int>     <dbl>     <dbl>     <dbl>
#>  1 AT1G01010.1 333554|PACid:16033839 0.106   0.254 0.420           73.2               347         474         75            8     52       369  77.8       1   430   430 100     100         1   466   466 1.37e-211       582      1500
#>  2 AT1G01020.1 470181|PACid:16064328 0.0402  0.104 0.388           91.1               224         246         22            0      0       231  93.9       1   246   246 100     100         1   246   246 2.27e-156       426      1094
#>  3 AT1G01030.1 470180|PACid:16054974 0.0150  0.126 0.118           93.3               335         359         20            2      4       338  94.2       1   359   359 100     100         1   355   355 4.01e-210       571      1471
#>  4 AT1G01040.1 333551|PACid:16057793 0.0135  0.116 0.116           93.4              1840        1969         58            7     71      1870  95         6  1910  1910  99.7    99.7       2  1963  1963 0              3544      9190
#>  5 AT1G01050.1 909874|PACid:16064489 0       0.175 0              100                 213         213          0            0      0       213 100         1   213   213 100     100         1   213   213 4.47e-158       427      1098
#>  6 AT1G01060.3 470177|PACid:16043374 0.0449  0.113 0.397           87.5               567         648         71            5     10       586  90.4       1   646   646 100     100         1   640   640 0              1028      2658
#>  7 AT1G01070.1 918864|PACid:16052578 0.0183  0.106 0.173           92.6               339         366         23            2      4       342  93.4       1   366   366 100     100         1   362   362 6.01e-227       614      1583
#>  8 AT1G01080.1 909871|PACid:16053217 0.0340  0.106 0.322           89.3               268         300         25            2      7       275  91.7       1   294   294 100     100         1   299   299 1.14e-181       494      1271
#>  9 AT1G01090.1 470171|PACid:16052860 0.00910 0.218 0.0417          96.8               420         434          8            3      6       425  97.9       1   429   429 100     100         1   433   433 5.41e-305       817      2110
#> 10 AT1G01110.2 333544|PACid:16034284 0.0325  0.122 0.266           87.7               463         528         65            0      0       475  90         1   528   528 100     100         1   528   528 1.18e-309       837      2162
#> 11 AT1G01120.1 918858|PACid:16049140 0.00307 0.133 0.0232          99.2               525         529          4            0      0       527  99.6       1   529   529 100     100         1   529   529 0              1029      2661
#> 12 AT1G01140.3 470161|PACid:16036015 0.00567 0.131 0.0432          98.5               446         453          6            1      1       448  98.9       1   452   452 100     100         1   453   453 0               865      2235
#> 13 AT1G01150.1 918855|PACid:16037307 0.13    0.203 0.641           72.6               207         285         68            3     10       234  82.1       5   289   346  82.4    82.4      16   290   294 5.31e-148       410      1055
#> 14 AT1G01160.1 918854|PACid:16044153 0.105   0.280 0.373           78.8               141         179         30            2      8       145  81         4   175   196  87.8    87.8       2   179   225 1.62e- 94       266       680
#> 15 AT1G01170.2 311317|PACid:16052302 0       0.306 0               85.6                83          97          0            1     14        83  85.6       2    84    84  98.8    98.8       4   100   100 8.03e- 57       161       408
#> 16 AT1G01180.1 909860|PACid:16056125 0.0297  0.176 0.168           92.6               287         310         20            1      3       296  95.5       1   307   479  64.1    64.1       1   310   316 4.57e-214       584      1506
#> 17 AT1G01190.1 311315|PACid:16059488 0.0287  0.162 0.177           94.2               502         533         30            1      1       518  97.2       1   533   536  99.4    99.4       6   537   540 0              1007      2603
#> 18 AT1G01200.1 470156|PACid:16041002 0.0190  0.168 0.114           95.8               228         238         10            0      0       233  97.9       1   238   238 100     100         1   238   238 6.96e-163       441      1135
#> 19 AT1G01210.1 311313|PACid:16057125 0.0207  0.154 0.134           95.3               102         107          5            0      0       105  98.1       1   107   107 100     100         1   107   107 2.04e- 80       222       566
#> 20 AT1G01220.1 470155|PACid:16047984 0.0157  0.153 0.102           96.5              1019        1056         37            0      0      1035  98         1  1056  1056 100     100         1  1056  1056 0              2016      5224
```

The [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
function can be used choosing the following options:

- `aligner` : `"diamond"` (default; DIAMOND2) or `"blast"` (BLASTP)
- `ortho_detection` : `RBH` (reciprocal best hit), `BH` (best hit) —
  works with both DIAMOND2 and BLAST depending on `aligner`
- `aa_aln_type` : `multiple` or `pairwise`
- `aa_aln_tool` : `clustalw`, `t_coffee`, `muscle`, `clustalo`, `mafft`,
  and `NW` (in case `aa_aln_type = "pairwise"`)
- `codon_aln_tool` : `pal2nal`
- `dnds_est.method` : `Li`, `Comeron`, `NG`, `LWL`, `LPB`, `MLWL`, `YN`,
  and `MYN`

Please see
[`?dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md) for
details.

In case your DIAMOND2 program or multiple alignment program cannot be
executed from the default execution `PATH` you can specify the
`aa_aln_path` or `aligner_path` arguments.

``` r

library(orthologr)
# using the `aa_aln_path` or `aligner_path` arguments
dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      aligner         = "diamond",
      ortho_detection = "RBH",
      aligner_path    = "here/path/to/diamond",
      aa_aln_type     = "multiple", 
      aa_aln_tool     = "clustalw", 
      aa_aln_path     = "here/path/to/clustalw",
      codon_aln_tool  = "pal2nal", 
      dnds_est.method = "Comeron", 
      comp_cores      = 1, 
      clean_folders   = TRUE)
```

## Advanced options

Additional arguments can be passed to
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md). This
allows you to use more advanced options of several interface programs.

To pass additional parameters to the interface programs, you can use the
`blast_params` and `aa_aln_params` arguments. The `aa_aln_params`
argument assumes that when you chose e.g. `aa_aln_tool = "mafft"` you
will pass the corresponding additional parameters in MAFFT notation.

``` r

library(orthologr)
# get dNdS estimated for orthologous genes between A. thaliana and A. lyrata
# using additional parameters:

# get a dNdS table using:
# 1) DIAMOND2 reciprocal best hit for orthology inference (RBH), with extra DIAMOND2 params
# 2) multiple amino acid alignments using MAFFT
# 3) pal2nal for codon alignments
# 4) Comeron (1995) for dNdS estimation
# 5) single core processing 'comp_cores = 1'
Ath_Aly_dnds <- dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
                      aligner         = "diamond",
                      ortho_detection = "RBH",
                      aligner_params  = "--matrix BLOSUM80",
                      aa_aln_type     = "multiple",
                      aa_aln_tool     = "mafft",
                      aa_aln_params   = "--maxiterate 1 --clustalout",
                      dnds_est.method = "Comeron",
                      comp_cores      = 1, 
                      clean_folders   = TRUE, 
                      quiet           = TRUE )

# filter for:
# 1) all dN values having an NA value are omitted
# 2) all dS values having an NA value are omitted
# 3) all dNdS values >= 0.1 are omitted
filter_dNdS(Ath_Aly_dnds, dnds.threshold = 0.1)
```


         query_id      subject_id       dN     dS    dNdS
    1 AT1G01050.1 909874|PACid:16 0.000000 0.1750 0.00000
    2 AT1G01090.1 470171|PACid:16 0.009843 0.2150 0.04579
    3 AT1G01120.1 918858|PACid:16 0.003072 0.1326 0.02317
    4 AT1G01140.3 470161|PACid:16 0.005672 0.1312 0.04324
    5 AT1G01170.2 311317|PACid:16 0.008750 0.2827 0.03095
    6 AT1G01220.1 470155|PACid:16 0.015210 0.1533 0.09919

Here `aligner_params` and `aa_aln_params` take a character string
specifying the parameters that shall be passed to DIAMOND2 (or BLAST
when `aligner = "blast"`) and MAFFT, respectively. The notation of these
parameters must follow the command line call of the stand alone versions
of DIAMOND2 and MAFFT: e.g. `aligner_params = "--matrix BLOSUM80"` and
`aa_aln_params = "--maxiterate 1 --clustalout"`.
