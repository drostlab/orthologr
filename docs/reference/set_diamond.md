# Create a DIAMONDable database with `diamond makedb`

This function reads a file storing a specific sequence type, such as
"cds", "protein", or "dna" in a standard sequence file format such as
"fasta", etc. and depending of the `makedb` parameter either creates a
diamond-able database, or returns the corresponding protein sequences as
data.table object for further DIAMOND2 searches.

## Usage

``` r
set_diamond(
  file,
  seq_type = "cds",
  format = "fasta",
  makedb = FALSE,
  delete_corrupt_cds = TRUE,
  path = NULL,
  makedb_type = "protein",
  comp_cores = 1,
  quiet = TRUE,
  ...
)
```

## Arguments

- file:

  a character string specifying the path to the file storing the
  sequences of interest.

- seq_type:

  a character string specifying the sequence type stored in the input
  file. Options are are: "cds", "protein", or "dna". In case of "cds",
  sequence are translated to protein sequences, in case of "dna", cds
  prediction is performed on the corresponding sequences which
  subsequently are translated to protein sequences. Default is
  `seq_type` = "cds".

- format:

  a character string specifying the file format used to store the
  genome, e.g. "fasta", "gbk".

- makedb:

  TRUE or FALSE whether a database should be created or not
  (`diamond makedb`).

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- path:

  a character string specifying the path to the DIAMOND2 program (in
  case you don't use the default path).

- makedb_type:

  a character string specifying the sequence type stored in the DIAMOND2
  database that is generated using 'diamond makedb'. Currently, the only
  option is "protein". Default is `makedb_type` = "protein".

- comp_cores:

  a numeric value specifying the number of cores to be used for
  multicore 'diamond makedb' computations.

- quiet:

  a logical value indicating whether `diamond makedb` should be run with
  the quiet mode. Default is `quiet` = `TRUE` (which adds `--quiet` to
  the diamond makedb run).

- ...:

  additional arguments that are used by the seqinr::read.fasta()
  function.

## Value

A list storing two elements. The first element \[\[1\]\] corresponds to
the data.table storing the gene ids in the first column and the
corresponding dna (cds) sequence in the second column and the aminoacid
sequence third column. The second list element \[\[2\]\] stores the name
of the protein database that was created by 'diamond makedb'.

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options

## See also

[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`set_blast`](https://drostlab.github.io/orthologr/reference/set_blast.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
 # running the set function to see an example output
 head(set_diamond(file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'))[[1]] , 2)
} # }
```
