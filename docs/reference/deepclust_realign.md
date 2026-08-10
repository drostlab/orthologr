# Realign sequences within `deepclust` clusters

Runs `diamond realign` on an existing set of
[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
clusters to obtain per-pair alignment statistics (approximate identity,
e-value, bitscore, and query/subject coverages) for every
`(member, representative)` pair inside each cluster. The result can be
passed directly to
[`deepclust_annotate`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
via its `realign` argument to enrich a cluster profile without
re-running `diamond deepclust`.

## Usage

``` r
deepclust_realign(
  input_file,
  clusters,
  seq_type = "protein",
  format = "fasta",
  delete_corrupt_cds = TRUE,
  path = NULL,
  comp_cores = 1,
  realign_params = NULL,
  save.output = NULL,
  quiet = TRUE
)
```

## Arguments

- input_file:

  a character string, or a character vector of file paths, specifying
  the input sequence file(s) that were used for the original
  [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  run. When multiple paths are provided they are merged into a single
  temporary FASTA before the DIAMOND2 database is built (same behaviour
  as
  [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)).
  Both plain and gzip-compressed (`.gz`) files are supported.

- clusters:

  either

  - a `character` string — path to a tab-separated `deepclust` output
    file whose first column is `representative_id` and second column is
    `member_id`, or

  - a `data.frame` / `tibble` with columns `representative_id` and
    `member_id` (i.e. the direct return value of
    [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
    or
    [`deepclust_annotate`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)).

- seq_type:

  a character string specifying the sequence type stored in the input
  file(s). Options are: `"cds"`, `"protein"`, or `"dna"`. In case of
  `"cds"`, sequences are translated to protein sequences before building
  the database. Default is `seq_type = "protein"`.

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

- realign_params:

  a character string of additional `diamond realign` parameters to
  append to the command. Default is `NULL`.

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
seven columns:

- `representative_id` — accession of the cluster representative (query
  in the realignment).

- `member_id` — accession of the cluster member (subject in the
  realignment).

- `approx_pident` — approximate percentage of identical matches.

- `evalue` — expect value.

- `bitscore` — bit score.

- `qcovhsp` — query coverage per HSP.

- `scovhsp` — subject (representative) coverage per HSP.

## Details

This function wraps `diamond realign`, which re-aligns every cluster
member against its assigned representative using the cluster mapping
produced by
[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md).
The output format is BLAST tabular (`--outfmt 6`) with the fields:
`qseqid sseqid approx_pident evalue bitscore qcovhsp scovhsp`.

The `--clusters` argument consumed by `diamond realign` expects a
two-column tab-separated file (`representative_id TAB member_id`). When
a tibble is supplied to `clusters`, it is written to a temporary file
automatically.

Columns in the returned tibble are named to match those used by
[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
(`representative_id`, `member_id`) so that the result joins cleanly onto
[`deepclust_annotate`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md)
output. Note that `diamond realign` uses the *representative* as the
alignment query (`qseqid`) and the *member* as the subject (`sseqid`);
the returned tibble is named accordingly.

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki

## See also

[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md),
[`deepclust_annotate`](https://drostlab.github.io/orthologr/reference/deepclust_annotate.md),
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. Run deepclust first
clusters <- deepclust(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  )
)

# 2. Realign each representative against its cluster members to get alignment stats
realign_result <- deepclust_realign(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  clusters = clusters
)

# 2. Alternatively, supply the path to the deepclust TSV directly
realign_result <- deepclust_realign(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  clusters = "/path/to/deepclust_output.tsv"
)

# 3. Incorporate into deepclust_annotate without re-running deepclust
profile <- deepclust_annotate(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  realign = realign_result
)
} # }
```
