# Compute dNdS values for two organisms

This function takes the CDS files of two organisms of interest
(query_file and subject_file) and computes the dNdS estimation values
for orthologous gene pairs between these organisms.

## Usage

``` r
dNdS(
  query_file,
  subject_file,
  aligner = "diamond",
  sensitivity_mode = "fast",
  aligner_path = NULL,
  seq_type = "cds",
  format = "fasta",
  ortho_detection = "RBH",
  delete_corrupt_cds = FALSE,
  store_locally = FALSE,
  cdd.path = NULL,
  aligner_params = NULL,
  eval = "1E-5",
  ortho_path = NULL,
  aa_aln_type = "pairwise",
  aa_aln_tool = "NW",
  aa_aln_path = NULL,
  aa_aln_params = NULL,
  codon_aln_tool = "pal2nal",
  kaks_calc_path = NULL,
  dnds_est.method = "Comeron",
  comp_cores = 1,
  quiet = TRUE,
  clean_folders = FALSE,
  print_citation = TRUE
)
```

## Arguments

- query_file:

  a character string specifying the path to the CDS file of interest
  (query organism).

- subject_file:

  a character string specifying the path to the CDS file of interest
  (subject organism).

- aligner:

  a character string specifying the sequence aligner. The options are
  `diamond` and `blast`. The option `diamond` uses DIAMOND2 which is up
  to 10 000X folds faster than BLASTP while retaining the sensitivity of
  BLASTP. Thus, the default is `aligner` = `diamond`.

- sensitivity_mode:

  specify the level of alignment sensitivity, when using DIAMOND2. The
  higher the sensitivity level, the more deep homologs can be found, but
  at the cost of reduced computational speed. - sensitivity_mode =
  "faster" : fastest alignment mode, but least sensitive (default).
  Designed for finding hits of \>70 - sensitivity_mode = "default" :
  Default mode. Designed for finding hits of \>70 - sensitivity_mode =
  "fast" : fast alignment mode, but least sensitive (default). Designed
  for finding hits of \>70 - sensitivity_mode = "mid-sensitive" : fast
  alignments between the fast mode and the sensitive mode in
  sensitivity. - sensitivity_mode = "sensitive" : fast alignments, but
  full sensitivity for hits \>40 - sensitivity_mode = "more-sensitive" :
  more sensitive than the sensitive mode. - sensitivity_mode =
  "very-sensitive" : sensitive alignment mode. - sensitivity_mode =
  "ultra-sensitive" : most sensitive alignment mode (sensitivity as high
  as BLASTP).

- aligner_path:

  a character string specifying the path to the DIAMOND or BLAST program
  (in case you don't use the default path).

- seq_type:

  a character string specifying the sequence type stored in the input
  file.Options are are:

  - `seq_type = "cds"` (Default): sequence are translated to protein
    sequences

  - `seq_type = "protein"`: orthology inference is performed using
    protein sequences directly.

- format:

  a character string specifying the file format of the sequence file,
  e.g. `format = "fasta"`, `format = "gbk"`. See
  [`read.cds`](https://drostlab.github.io/orthologr/reference/read.cds.md),
  [`read.genome`](https://drostlab.github.io/orthologr/reference/read.genome.md),
  [`read.proteome`](https://drostlab.github.io/orthologr/reference/read.proteome.md)
  for more details.

- ortho_detection:

  a character string specifying the orthology inference method that
  shall be performed to detect orthologous genes. Options are:

  - `ortho_detection ="BH"`: DIAMOND or BLAST unidirectional best hit.

  - `ortho_detection = "RBH"`: DIAMOND or BLAST reciprocal/bidirectional
    best hit (Default).

  - `ortho_detection = "Orthofinder2"`: single copy core orthologs
    between multiple species proteome comparisons.

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- store_locally:

  a logical value indicating whether or not alignment files shall be
  stored locally rather than in
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- cdd.path:

  path to the cdd database folder (specify when using `ortho_detection`
  = `"DELTA"`).

- aligner_params:

  a character string specifying additional parameters that shall be
  passed to DIAMOND or BLAST. Default is `aligner_params` = `NULL`.

- eval:

  a numeric value specifying the E-Value cutoff for DIAMOND or BLAST hit
  detection.

- ortho_path:

  a character string specifying the path to the orthology inference
  program such as `Orthofinder2`, etc. (in case you don't use the
  default path).

- aa_aln_type:

  a character string specifying the amino acid alignment type:

  - `aa_aln_type` = "multiple" (Default)

  - `aa_aln_type` = "pairwise"

  .

- aa_aln_tool:

  a character string specifying the program that should be used e.g.
  "clustalw".

- aa_aln_path:

  a character string specifying the path to the multiple alignment
  program (in case you don't use the default path).

- aa_aln_params:

  a character string specifying additional parameters that shall be
  passed to the selected alignment tool. Default is `aa_aln_params` =
  `NULL` (no addintional parameters are passed to the selected alignment
  tool).

- codon_aln_tool:

  a character string specifying the codon alignment tool that shall be
  used. Default is `codon_aln_tool` = `"pal2nal"`. Right now only
  "pal2nal" can be selected as codon alignment tool.

- kaks_calc_path:

  a character string specifying the execution path to KaKs_Calculator.
  Default is `kaks_calc_path` = `NULL` (meaning that KaKs_Calculator is
  stored and executable in your default `PATH`).

- dnds_est.method:

  the dNdS estimation method that shall be used. Options are:

  - `dnds_est.method = "Comeron"` (Default): Comeron's method (1995)

  - `dnds_est.method = "Li"`: Li's method (1993)

  - `dnds_est.method = "NG"`: Nei, M. and Gojobori, T. (1986)

  - `dnds_est.method = "LWL"`: Li, W.H., et al. (1985)

  - `dnds_est.method = "LPB"`: Li, W.H. (1993) and Pamilo, P. and
    Bianchi, N.O. (1993)

  - `dnds_est.method = "MLWL"`: (Modified LWL), MLPB (Modified LPB):
    Tzeng, Y.H., et al. (2004)

  - `dnds_est.method = "YN"`: Yang, Z. and Nielsen, R. (2000)

  - `dnds_est.method ="MYN"` (Modified YN): Zhang, Z., et al. (2006)

- comp_cores:

  a numeric value specifying the number of cores that shall be used to
  perform parallel computations on a multicore machine.

- quiet:

  a logical value specifying whether the output of the corresponding
  alignment tool shall be printed out to the console. Default is `quiet`
  = `FALSE`.

- clean_folders:

  a boolean value spefiying whether all internall folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

- print_citation:

  a logical value indicating whether or not the citation message shall
  be printed.

## Value

A data.table storing the dNdS values of the correspnding genes.

## Details

The dN/dS ratio quantifies the mode and strength of selection acting on
a pair of orthologous genes. This selection pressure can be quantified
by comparing synonymous substitution rates (dS) that are assumed to be
neutral with nonsynonymous substitution rates (dN), which are exposed to
selection as they change the amino acid composition of a protein (Mugal
et al., 2013 http://mbe.oxfordjournals.org/content/31/1/212).

The orthologr package provides the `dNdS` function to perform dNdS
estimation on pairs of orthologous genes. This function takes the CDS
files of two organisms of interest (`query_file` and `subject_file`) and
computes the dNdS estimation values for orthologous gene pairs between
these organisms.

The following pipieline resembles the dNdS estimation process:

- 1\) Orthology Inference: e.g. DIAMOND or BLAST reciprocal best hit
  (RBH)

- 2\) Pairwise sequence alignment: e.g. clustalw for pairwise amino acid
  sequence alignments

- 3\) Codon Alignment: e.g. pal2nal program

- 4\) dNdS estimation: e.g. Yang, Z. and Nielsen, R. (2000)
  <http://mbe.oxfordjournals.org/content/17/1/32.short>

Note: it is assumed that when using `dNdS()` all corresponding multiple
sequence alignment programs you want to use are already installed on
your machine and are executable via either the default execution `PATH`
or you specifically define the location of the executable file via the
`aa_aln_path` or `aligner_path` argument that can be passed to `dNdS()`.

The `dNdS()` function can be used choosing the folllowing options:

- `ortho_detection` :

  - `"RBH"` (DIAMOND or BLAST best reciprocal hit)

  - `"BH"` (DIAMOND or BLAST best hit)

  - `"Orthofinder2"`

- `aa_aln_type` :

  - `"multiple"`

  - `"pairwise"`

- `aa_aln_tool` :

  - `"clustalw"`

  - `"t_coffee"`

  - `"muscle"`

  - `"clustalo"`

  - `"mafft"`

  - `"NW"` (in case `aa_aln_type = "pairwise"`)

- `codon_aln_tool` :

  - `"pal2nal"`

- `dnds_est.method` :

  - "Li" : Li's method (1993)

  - "Comeron" : Comeron's method (1995)

  - "NG": Nei, M. and Gojobori, T. (1986)

  - "LWL": Li, W.H., et al. (1985)

  - "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)

  - "MLWL" (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al.
    (2004)

  - "YN": Yang, Z. and Nielsen, R. (2000)

  - "MYN" (Modified YN): Zhang, Z., et al. (2006)

## References

seqinr: <https://cran.r-project.org/package=seqinr>

Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator:
Calculating Ka and Ks through model selection and model averaging.
Genomics Proteomics Bioinformatics 2006 , 4:259-263.

KaKs_Calculator: <https://code.google.com/archive/p/kaks-calculator>
\[GNU GPL-3 license\]

Paradis, E. (2012) Analysis of Phylogenetics and Evolution with R
(Second Edition). New York: Springer.

Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of
phylogenetics and evolution in R language. Bioinformatics, 20, 289-290.

More information on ape can be found at
<https://emmanuelparadis.github.io/> or
<https://github.com/emmanuelparadis/ape>.

Pages H, Aboyoun P, Gentleman R and DebRoy S. Biostrings: String objects
representing biological sequences, and matching algorithms. R package
version 2.32.1.

## See also

[`divergence_stratigraphy`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md),
[`orthologs`](https://drostlab.github.io/orthologr/reference/orthologs.md),
[`substitutionrate`](https://drostlab.github.io/orthologr/reference/substitutionrate.md),
[`multi_aln`](https://drostlab.github.io/orthologr/reference/multi_aln.md),
[`codon_aln`](https://drostlab.github.io/orthologr/reference/codon_aln.md),
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`blast_best`](https://drostlab.github.io/orthologr/reference/blast_best.md),
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`read.cds`](https://drostlab.github.io/orthologr/reference/read.cds.md)

## Author

Hajk-Georg Drost and Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{

# get a dNdS table using:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single core processing 'comp_cores = 1'
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1 )


# running dNdS using the 'aa_aln_path' argument to specify the path to
# the corresponding alignment tool
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     aa_aln_path     = "/path/to/clustalw/",
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1 )

# The same result can be obtained using multicore processing using: comp_cores = 2 or 3 or more ...
dNdS(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
     subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     aa_aln_path     = "/path/to/clustalw/",
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1 )

} # }
```
