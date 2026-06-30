# Cluster sequences using `diamond deepclust`

This function clusters a set of protein (or CDS) sequences using the
DIAMOND2 `deepclust` algorithm. It prepares a DIAMOND2 database via
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
runs `diamond deepclust`, and returns a two-column tibble mapping each
representative sequence to its cluster members. The result is also
written to a TSV file. Note, we recommend users to run
[`deepclust_annotate()`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
for more information about the cluster members and their representative
sequences.

## Usage

``` r
deepclust(
  input_file,
  seq_type = "protein",
  format = "fasta",
  delete_corrupt_cds = TRUE,
  path = NULL,
  comp_cores = 1,
  mutual_cover = NULL,
  approx_id = NULL,
  eval = NULL,
  deepclust_params = NULL,
  save.output = NULL,
  quiet = TRUE
)
```

## Arguments

- input_file:

  a character string, or a character vector of file paths, specifying
  the input sequence file(s). When multiple paths are provided the
  sequences are merged into a single temporary FASTA before clustering,
  allowing proteomes from different organisms to be clustered together.
  Both plain and gzip-compressed (`.gz`) files are supported and may be
  mixed freely. Supported sequence types are controlled by `seq_type`.

- seq_type:

  a character string specifying the sequence type stored in the input
  file(s). Options are: `"cds"`, `"protein"`, or `"dna"`. In case of
  `"cds"`, sequences are translated to protein sequences before
  clustering. Default is `seq_type = "protein"`.

- format:

  a character string specifying the file format of the sequence file.
  Default is `format = "fasta"`.

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input file. Only relevant when
  `seq_type = "cds"`. Default is `delete_corrupt_cds = TRUE`.

- path:

  a character string specifying the path to the DIAMOND2 executable (in
  case it is not on the system `PATH`).

- comp_cores:

  a numeric value specifying the number of CPU threads to use. Passed to
  `--threads`. Default is `comp_cores = 1`.

- mutual_cover:

  a numeric value (percentage) specifying the minimum mutual coverage of
  the cluster member and representative sequence. Passed to
  `--mutual-cover`. Default is `NULL` (DIAMOND2 default applies).

- approx_id:

  a numeric value (percentage) specifying the minimum approximate
  identity to report an alignment or to cluster sequences. Passed to
  `--approx-id`. Default is `NULL` (DIAMOND2 default applies).

- eval:

  a numeric value specifying the maximum E-value to report alignments.
  Passed to `--evalue`. Default is `NULL` (DIAMOND2 default of 0.001
  applies).

- deepclust_params:

  a character string of additional DIAMOND2 deepclust parameters to
  append to the `deepclust` call. Default is `NULL`.

- save.output:

  a character string specifying the path where the output TSV file
  should be saved. E.g. `save.output = getwd()` to save in the current
  working directory. Default is `NULL` (file is only written to a
  temporary directory).

- quiet:

  a logical value indicating whether DIAMOND2 should run in quiet mode.
  Default is `quiet = TRUE`.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
two columns:

- `representative_id` — accession of the cluster representative
  sequence.

- `member_id` — accession of the cluster member sequence.

## Details

This function wraps `diamond deepclust`, the graph-based sequence
clustering algorithm available in DIAMOND2 v2.1.0 and later. The output
is a two-column tab-separated file where the first column contains the
representative (centroid) sequence accession and the second column
contains the cluster member accession. A sequence that is its own
representative will appear with itself in both columns.

When multiple files are supplied via `input_file`, they are first merged
into a single temporary FASTA file (using
[`seqinr::read.fasta`](https://rdrr.io/pkg/seqinr/man/read.fasta.html)
and
[`seqinr::write.fasta`](https://rdrr.io/pkg/seqinr/man/write.fasta.html))
before the DIAMOND2 database is built. This allows cross-species or
multi-proteome clustering in a single run. The merged file is written to
the session temporary directory and is not retained after the R session
ends unless `save.output` is specified.

The function uses
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md)
internally to handle CDS translation and DIAMOND2 database creation.

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki

## See also

[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md),
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`deepclust_annotate`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md),
[`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
# Cluster a single proteome (protein sequences, the default)
deepclust(
  input_file = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr')
)

# Cluster two proteomes together in a single run
deepclust(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  )
)

# Apply identity and coverage thresholds across multiple proteomes
deepclust(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  approx_id    = 50,
  mutual_cover = 80,
  eval         = 1e-5,
  comp_cores   = 4
)

# Save the cluster TSV to the current working directory
deepclust(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  save.output = getwd()
)

# Cluster CDS sequences (translated to protein internally)
deepclust(
  input_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
  seq_type   = "cds"
)
} # }
```
