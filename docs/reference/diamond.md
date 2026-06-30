# Perform a DIAMOND2 search

This function performs a DIAMOND2 search of a given set of sequences
against a given database.

## Usage

``` r
diamond(
  query_file,
  subject_file,
  seq_type = "cds",
  format = "fasta",
  diamond_algorithm = "blastp",
  sensitivity_mode = "fast",
  eval = "1E-5",
  max.target.seqs = 10000,
  delete_corrupt_cds = TRUE,
  path = NULL,
  comp_cores = 1,
  diamond_params = NULL,
  clean_folders = FALSE,
  save.output = NULL,
  quiet = TRUE,
  database_maker = "diamond"
)
```

## Arguments

- query_file:

  a character string specifying the path to the CDS file of interest
  (query organism).

- subject_file:

  a character string specifying the path to the CDS file of interest
  (subject organism).

- seq_type:

  a character string specifying the sequence type stored in the input
  file. Options are are: "cds", "protein", or "dna". In case of "cds",
  sequence are translated to protein sequences, in case of "dna", cds
  prediction is performed on the corresponding sequences which
  subsequently are translated to protein sequences. Default is
  `seq_type` = "cds".

- format:

  a character string specifying the file format of the sequence file,
  e.g. `format` = `"fasta"`. Default is `format` = `"fasta"`.

- diamond_algorithm:

  a character string specifying the DIAMOND2 algorithm that shall be
  used, option is currently limited to: `diamond_algorithm` = `"blastp"`

- sensitivity_mode:

  specify the level of alignment sensitivity. The higher the sensitivity
  level, the more deep homologs can be found, but at the cost of reduced
  computational speed. - sensitivity_mode = "faster" : fastest alignment
  mode, but least sensitive (default). Designed for finding hits of
  \>70 - sensitivity_mode = "default" : Default mode. Designed for
  finding hits of \>70 - sensitivity_mode = "fast" : fast alignment
  mode, but least sensitive (default). Designed for finding hits of
  \>70 - sensitivity_mode = "mid-sensitive" : fast alignments between
  the fast mode and the sensitive mode in sensitivity. -
  sensitivity_mode = "sensitive" : fast alignments, but full sensitivity
  for hits \>40 - sensitivity_mode = "more-sensitive" : more sensitive
  than the sensitive mode. - sensitivity_mode = "very-sensitive" :
  sensitive alignment mode. - sensitivity_mode = "ultra-sensitive" :
  most sensitive alignment mode (sensitivity as high as BLASTP).

- eval:

  a numeric value specifying the E-Value cutoff for DIAMOND2 hit
  detection.

- max.target.seqs:

  a numeric value specifying the number of aligned sequences to keep.
  Please be aware that `max.target.seqs` selects best hits based on the
  database entry and not by the best e-value. See details here:
  https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166
  .

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- path:

  a character string specifying the path to the DIAMOND2 program (in
  case you don't use the default path).

- comp_cores:

  a numeric value specifying the number of cores that shall be used to
  run DIAMOND2 searches.

- diamond_params:

  a character string listing the input parameters that shall be passed
  to the executing DIAMOND2 program. Default is `NULL`, implicating that
  a set of default parameters is used when running DIAMOND2.

- clean_folders:

  a boolean value specifying whether all internal folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

- save.output:

  a path to the location were the DIAMOND2 output shall be stored. E.g.
  `save.output` = [`getwd()`](https://rdrr.io/r/base/getwd.html) to
  store it in the current working directory, or `save.output` =
  `file.path(put,your,path,here)`.

- quiet:

  a logical value indicating whether DIAMOND2 should be run with the
  quiet mode. Default is `quiet` = `TRUE` (which adds `--quiet` to the
  diamond run).

- database_maker:

  a character string specifying whether the database should be made
  using diamond or blast. Default is `database_maker` = `diamond`.

## Value

A data.table storing the DIAMOND2 hit table returned by DIAMOND2. The
format is the same as with BLAST.

## Details

This function provides a fast communication between R and DIAMOND2. It
is mainly used as internal functions such as
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md)
and
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md)
but can also be used to perform simple DIAMOND2 computations. This
function gives the same output as
[`blast`](https://drostlab.github.io/orthologr/reference/blast.md) while
being up to 10 000X faster in larger databases.

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options

## See also

[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
[`blast`](https://drostlab.github.io/orthologr/reference/blast.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
# performing a DIAMOND2 search using diamond blastp (default)
diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))

# performing a DIAMOND2 search using diamond blastp (default) 
# using amino acid sequences as input file
diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
        subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
        seq_type     = "protein")


# save the DIAMOND2 output table in your current working directory
diamond(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
        subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
        seq_type     = "protein",
        save.output  = getwd())

# in case you are working with a multicore machine, you can also run parallel
# DIAMOND2 computations using the comp_cores parameter: here with 2 cores
diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
        subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
        comp_cores   = 2)


 # running diamond using additional parameters
 diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
         diamond_params = "--max-target-seqs 1")


# running diamond using additional parameters and an external diamond path
 diamond(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
         subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
         diamond_params = "--max-target-seqs 1", path = "path/to/diamond/")
} # }
```
