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
```

    Running diamond version 2.2.1 ...
    sensitivity mode: fast
    creating a diamond database
    Starting DIAMOND2 search ...
    DIAMOND2 search completed in 0.02 sec.
    Running diamond version 2.2.1 ...
    sensitivity mode: fast
    creating a diamond database
    Starting DIAMOND2 search ...
    DIAMOND2 search completed in 0.02 sec.
    # A tibble: 20 × 21
    # Groups:   query_id [20]
       query_id    subject_id        perc_identity num_ident_matches alig_length
       <chr>       <chr>                     <dbl>             <int>       <int>
     1 AT1G01010.1 333554|PACid:160…          73.2               347         474
     2 AT1G01020.1 470181|PACid:160…          91.1               224         246
     3 AT1G01030.1 470180|PACid:160…          93.3               335         359
     4 AT1G01040.1 333551|PACid:160…          93.4              1840        1969
     5 AT1G01050.1 909874|PACid:160…         100                 213         213
     6 AT1G01060.3 470177|PACid:160…          87.5               567         648
     7 AT1G01070.1 918864|PACid:160…          92.6               339         366
     8 AT1G01080.1 909871|PACid:160…          89.3               268         300
     9 AT1G01090.1 470171|PACid:160…          96.8               420         434
    10 AT1G01110.2 333544|PACid:160…          87.7               463         528
    11 AT1G01120.1 918858|PACid:160…          99.2               525         529
    12 AT1G01140.3 470161|PACid:160…          98.5               446         453
    13 AT1G01150.1 918855|PACid:160…          72.6               207         285
    14 AT1G01160.1 918854|PACid:160…          78.8               141         179
    15 AT1G01170.2 311317|PACid:160…          85.6                83          97
    16 AT1G01180.1 909860|PACid:160…          92.6               287         310
    17 AT1G01190.1 311315|PACid:160…          94.2               502         533
    18 AT1G01200.1 470156|PACid:160…          95.8               228         238
    19 AT1G01210.1 311313|PACid:160…          95.3               102         107
    20 AT1G01220.1 470155|PACid:160…          96.5              1019        1056
    # ℹ 16 more variables: mismatches <int>, gap_openings <int>, n_gaps <int>,
    #   pos_match <int>, ppos <dbl>, q_start <int>, q_end <int>, q_len <int>,
    #   qcov <dbl>, qcovhsp <dbl>, s_start <int>, s_end <int>, s_len <int>,
    #   evalue <dbl>, bit_score <dbl>, score_raw <dbl>

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
```

    # A tibble: 20 x 21
    # Groups:   query_id [20]
       query_id subject_id perc_identity alig_length
       <chr>    <chr>              <dbl>       <dbl>
     1 AT1G010... 333554|PA...          74.0         469
     2 AT1G010... 470181|PA...          91.1         246
     3 AT1G010... 470180|PA...          95.5         359
     4 AT1G010... 333551|PA...          92.0        1970
     5 AT1G010... 909874|PA...         100           213
     6 AT1G010... 470177|PA...          89.5         648
     7 AT1G010... 918864|PA...          95.1         366
     8 AT1G010... 909871|PA...          90.3         300
     9 AT1G010... 470171|PA...          96.8         434
    10 AT1G011... 333544|PA...          93.6         528
    11 AT1G011... 918858|PA...          99.2         529
    12 AT1G011... 470161|PA...          98.4         453
    13 AT1G011... 918855|PA...          72.6         285
    14 AT1G011... 918854|PA...          84.9         179
    15 AT1G011... 311317|PA...          85.6          97
    16 AT1G011... 909860|PA...          92.6         310
    17 AT1G011... 311315|PA...          94.2         533
    18 AT1G012... 470156|PA...          95.8         238
    19 AT1G012... 311313|PA...          95.3         107
    20 AT1G012... 470155|PA...          96.7        1056
    # ... with 17 more variables: mismatches <dbl>,
    #   gap_openings <dbl>, n_gaps <dbl>, pos_match <dbl>, ppos <dbl>,
    #   q_start <dbl>, q_end <dbl>, q_len <dbl>, qcov <dbl>, qcovhsp <dbl>,
    #   s_start <dbl>, s_end <dbl>, s_len <dbl>, evalue <dbl>,
    #   bit_score <dbl>, score_raw <dbl>

In this case 2 cores are being used to perform parallel processing, the
`clean_folders` argument specifies that all files returned by the
corresponding orthology inference method are removed after analyses.

## Program specific use of the orthologs() function

In this section small examples will illustrate the use of the
[`orthologs()`](https://drostlab.github.io/orthologr/reference/orthologs.md)
function for each orthology inference program.

### DIAMOND2 best hit (recommended)

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
```

    # A tibble: 20 × 21
    # Groups:   query_id [20]
       query_id    subject_id  perc_identity num_ident_matches alig_length
       <chr>       <chr>               <dbl>             <int>       <int>
     1 AT1G01010.1 333554|PAC…          73.2               347         474
     2 AT1G01020.1 470181|PAC…          91.1               224         246
     3 AT1G01030.1 470180|PAC…          93.3               335         359
     4 AT1G01040.1 333551|PAC…          93.4              1840        1969
     5 AT1G01050.1 909874|PAC…         100                 213         213
     6 AT1G01060.3 470177|PAC…          87.5               567         648
     7 AT1G01070.1 918864|PAC…          92.6               339         366
     8 AT1G01080.1 909871|PAC…          89.3               268         300
     9 AT1G01090.1 470171|PAC…          96.8               420         434
    10 AT1G01110.2 333544|PAC…          87.7               463         528
    11 AT1G01120.1 918858|PAC…          99.2               525         529
    12 AT1G01140.3 470161|PAC…          98.5               446         453
    13 AT1G01150.1 918855|PAC…          72.6               207         285
    14 AT1G01160.1 918854|PAC…          78.8               141         179
    15 AT1G01170.2 311317|PAC…          85.6                83          97
    16 AT1G01180.1 909860|PAC…          92.6               287         310
    17 AT1G01190.1 311315|PAC…          94.2               502         533
    18 AT1G01200.1 470156|PAC…          95.8               228         238
    19 AT1G01210.1 311313|PAC…          95.3               102         107
    20 AT1G01220.1 470155|PAC…          96.5              1019        1056
    # ℹ 16 more variables: mismatches <int>, gap_openings <int>,
    #   n_gaps <int>, pos_match <int>, ppos <dbl>, q_start <int>,
    #   q_end <int>, q_len <int>, qcov <dbl>, qcovhsp <dbl>,
    #   s_start <int>, s_end <int>, s_len <int>, evalue <dbl>,
    #   bit_score <dbl>, score_raw <dbl>

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
```

    # A tibble: 20 × 21
    # Groups:   query_id [20]
       query_id    subject_id  perc_identity num_ident_matches alig_length
       <chr>       <chr>               <dbl>             <int>       <int>
     1 AT1G01010.1 333554|PAC…          73.2               347         474
     2 AT1G01020.1 470181|PAC…          91.1               224         246
     3 AT1G01030.1 470180|PAC…          93.3               335         359
     4 AT1G01040.1 333551|PAC…          93.4              1840        1969
     5 AT1G01050.1 909874|PAC…         100                 213         213
     6 AT1G01060.3 470177|PAC…          87.5               567         648
     7 AT1G01070.1 918864|PAC…          92.6               339         366
     8 AT1G01080.1 909871|PAC…          89.3               268         300
     9 AT1G01090.1 470171|PAC…          96.8               420         434
    10 AT1G01110.2 333544|PAC…          87.7               463         528
    11 AT1G01120.1 918858|PAC…          99.2               525         529
    12 AT1G01140.3 470161|PAC…          98.5               446         453
    13 AT1G01150.1 918855|PAC…          72.6               207         285
    14 AT1G01160.1 918854|PAC…          78.8               141         179
    15 AT1G01170.2 311317|PAC…          85.6                83          97
    16 AT1G01180.1 909860|PAC…          92.6               287         310
    17 AT1G01190.1 311315|PAC…          94.2               502         533
    18 AT1G01200.1 470156|PAC…          95.8               228         238
    19 AT1G01210.1 311313|PAC…          95.3               102         107
    20 AT1G01220.1 470155|PAC…          96.5              1019        1056
    # ℹ 16 more variables: mismatches <int>, gap_openings <int>,
    #   n_gaps <int>, pos_match <int>, ppos <dbl>, q_start <int>,
    #   q_end <int>, q_len <int>, qcov <dbl>, qcovhsp <dbl>,
    #   s_start <int>, s_end <int>, s_len <int>, evalue <dbl>,
    #   bit_score <dbl>, score_raw <dbl>

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

# perform orthology inference using DIAMOND2 reciprocal best hit
RBH <- orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_RBH",
          clean_folders   = FALSE,
          comp_cores      = 1)

glimpse(RBH)
```

    Rows: 20
    Columns: 21
    Groups: query_id [20]
    $ query_id          <chr> "AT1G01010.1", "AT1G01020.1", "AT1G01030.1…
    $ subject_id        <chr> "333554|PACid:16033839", "470181|PACid:160…
    $ perc_identity     <dbl> 73.2, 91.1, 93.3, 93.4, 100.0, 87.5, 92.6,…
    $ num_ident_matches <int> 347, 224, 335, 1840, 213, 567, 339, 268, 4…
    $ alig_length       <int> 474, 246, 359, 1969, 213, 648, 366, 300, 4…
    $ mismatches        <int> 75, 22, 20, 58, 0, 71, 23, 25, 8, 65, 4, 6…
    $ gap_openings      <int> 8, 0, 2, 7, 0, 5, 2, 2, 3, 0, 0, 1, 3, 2, …
    $ n_gaps            <int> 52, 0, 4, 71, 0, 10, 4, 7, 6, 0, 0, 1, 10,…
    $ pos_match         <int> 369, 231, 338, 1870, 213, 586, 342, 275, 4…
    $ ppos              <dbl> 77.8, 93.9, 94.2, 95.0, 100.0, 90.4, 93.4,…
    $ q_start           <int> 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 5, 4, …
    $ q_end             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 4…
    $ q_len             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 4…
    $ qcov              <dbl> 100.0, 100.0, 100.0, 99.7, 100.0, 100.0, 1…
    $ qcovhsp           <dbl> 100.0, 100.0, 100.0, 99.7, 100.0, 100.0, 1…
    $ s_start           <int> 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 16, 2,…
    $ s_end             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 4…
    $ s_len             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 4…
    $ evalue            <dbl> 1.37e-211, 2.27e-156, 4.01e-210, 0.00e+00,…
    $ bit_score         <dbl> 582, 426, 571, 3544, 427, 1028, 614, 494, …
    $ score_raw         <dbl> 1500, 1094, 1471, 9190, 1098, 2658, 1583, …

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
```

    # A tibble: 20 × 21
    # Groups:   query_id [20]
       query_id    subject_id  perc_identity num_ident_matches alig_length
       <chr>       <chr>               <dbl>             <int>       <int>
     1 AT1G01010.1 333554|PAC…          74.0               347         469
     2 AT1G01020.1 470181|PAC…          91.1               224         246
     3 AT1G01030.1 470180|PAC…          95.5               343         359
     4 AT1G01040.1 333551|PAC…          92.0              1812        1970
     5 AT1G01050.1 909874|PAC…         100                 213         213
     6 AT1G01060.3 470177|PAC…          89.5               580         648
     7 AT1G01070.1 918864|PAC…          95.1               348         366
     8 AT1G01080.1 909871|PAC…          90.3               271         300
     9 AT1G01090.1 470171|PAC…          96.8               420         434
    10 AT1G01110.2 333544|PAC…          93.6               494         528
    11 AT1G01120.1 918858|PAC…          99.2               525         529
    12 AT1G01140.3 470161|PAC…          98.5               446         453
    13 AT1G01150.1 918855|PAC…          72.6               207         285
    14 AT1G01160.1 918854|PAC…          84.9               152         179
    15 AT1G01170.2 311317|PAC…          85.6                83          97
    16 AT1G01180.1 909860|PAC…          92.6               287         310
    17 AT1G01190.1 311315|PAC…          94.2               502         533
    18 AT1G01200.1 470156|PAC…          95.8               228         238
    19 AT1G01210.1 311313|PAC…          95.3               102         107
    20 AT1G01220.1 470155|PAC…          96.7              1021        1056
    # ℹ 16 more variables: mismatches <int>, gap_openings <int>,
    #   n_gaps <int>, pos_match <int>, ppos <dbl>, q_start <int>,
    #   q_end <int>, q_len <int>, qcov <dbl>, qcovhsp <dbl>,
    #   s_start <int>, s_end <int>, s_len <int>, evalue <dbl>,
    #   bit_score <dbl>, score_raw <dbl>

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
```

    # A tibble: 20 × 21
    # Groups:   query_id [20]
       query_id    subject_id  perc_identity num_ident_matches alig_length
       <chr>       <chr>               <dbl>             <int>       <int>
     1 AT1G01010.1 333554|PAC…          74.0               347         469
     2 AT1G01020.1 470181|PAC…          91.1               224         246
     3 AT1G01030.1 470180|PAC…          95.5               343         359
     4 AT1G01040.1 333551|PAC…          92.0              1812        1970
     5 AT1G01050.1 909874|PAC…         100                 213         213
     6 AT1G01060.3 470177|PAC…          89.5               580         648
     7 AT1G01070.1 918864|PAC…          95.1               348         366
     8 AT1G01080.1 909871|PAC…          90.3               271         300
     9 AT1G01090.1 470171|PAC…          96.8               420         434
    10 AT1G01110.2 333544|PAC…          93.6               494         528
    11 AT1G01120.1 918858|PAC…          99.2               525         529
    12 AT1G01140.3 470161|PAC…          98.5               446         453
    13 AT1G01150.1 918855|PAC…          72.6               207         285
    14 AT1G01160.1 918854|PAC…          84.9               152         179
    15 AT1G01170.2 311317|PAC…          85.6                83          97
    16 AT1G01180.1 909860|PAC…          92.6               287         310
    17 AT1G01190.1 311315|PAC…          94.2               502         533
    18 AT1G01200.1 470156|PAC…          95.8               228         238
    19 AT1G01210.1 311313|PAC…          95.3               102         107
    20 AT1G01220.1 470155|PAC…          96.7              1021        1056
    # ℹ 16 more variables: mismatches <int>, gap_openings <int>,
    #   n_gaps <int>, pos_match <int>, ppos <dbl>, q_start <int>,
    #   q_end <int>, q_len <int>, qcov <dbl>, qcovhsp <dbl>,
    #   s_start <int>, s_end <int>, s_len <int>, evalue <dbl>,
    #   bit_score <dbl>, score_raw <dbl>

In case you would like to store the corresponding `hit tables` returned
by BLAST for subsequent analyses, you can specify the
`clean_folders = FALSE` argument. The corresponding BLAST hit table can
then be found in `file.path(tempdir(),"_blast_db")`.

A detailed overview of further analyses that can be done with the
corresponding BLAST output can be found in the [BLAST
vignette](https://drostlab.github.io/orthologr/articles/blast.html).
