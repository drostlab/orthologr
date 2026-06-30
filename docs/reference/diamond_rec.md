# Perform a DIAMOND2 reciprocal best hit (RBH) search

This function performs a DIAMOND2 search (reciprocal best hit) of a
given set of protein sequences against a second set of protein sequences
and vice versa.

## Usage

``` r
diamond_rec(
  query_file,
  subject_file,
  seq_type = "cds",
  format = "fasta",
  diamond_algorithm = "blastp",
  sensitivity_mode = "fast",
  delete_corrupt_cds = TRUE,
  eval = "1E-5",
  max.target.seqs = 10000,
  path = NULL,
  comp_cores = 1,
  diamond_params = NULL,
  clean_folders = FALSE,
  save.output = NULL
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

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- eval:

  a numeric value specifying the E-Value cutoff for DIAMOND2 hit
  detection.

- max.target.seqs:

  a numeric value specifying the number of aligned sequences to keep.
  Please be aware that `max.target.seqs` selects best hits based on the
  database entry and not by the best e-value. See details here:
  https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166
  .

- path:

  a character string specifying the path to the DIAMOND2 program (in
  case you don't use the default path).

- comp_cores:

  a numeric value specifying the number of cores to be used for
  multicore DIAMOND2 computations.

- diamond_params:

  a character string listing the input paramters that shall be passed to
  the executing DIAMOND2 program. Default is `NULL`, implicating that a
  set of default parameters is used when running DIAMOND2.

- clean_folders:

  a boolean value spefiying whether all internall folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

- save.output:

  a path to the location were the DIAMOND2 output shall be stored. E.g.
  `save.output` = [`getwd()`](https://rdrr.io/r/base/getwd.html) to
  store it in the current working directory, or `save.output` =
  `file.path(put,your,path,here)`.

## Value

A data.table as returned by the
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md)
function, storing the geneids of orthologous genes (reciprocal best hit)
in the first column and the amino acid sequences in the second column.

## Details

Given a set of protein sequences A and a different set of protein
sequences B, first a best hit diamond search is being performed from A
to B: diamond(A,B) and afterwards a best hit diamond search is being
performed from B to A: diamond(B,A). Only protein sequences that were
found to be best hits in both directions are retained and returned.

This function gives the same output as
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md)
while being much much faster.

This function can be used to perform orthology inference using DIAMOND2
best reciprocal hit methodology.

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options

## See also

[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
# performing gene orthology inference using the reciprocal best hit (RBH) method
diamond_rec(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'))
          
          
          
# performing gene orthology inference using the reciprocal best hit (RBH) method
# starting with protein sequences
diamond_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
            seq_type     = "protein")



# save the DIAMOND2 output file to the current working directory
diamond_rec(query_file   = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
            subject_file = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
            seq_type     = "protein",
            save.output  = getwd())



# use multicore processing
diamond_rec(query_file    = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'), 
            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
            comp_cores   = 2)
           
           
           
# performing gene orthology inference using the reciprocal best hit (RBH) method
# and external path to diamond
diamond_rec(query_file   = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
            subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
            path         = "path/to/diamond")
          
          
          
} # }
```
