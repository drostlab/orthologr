# Annotate `deepclust` output with source file metadata

Runs
[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
across one or more sequence files, then annotates each cluster member
with the source file it originated from. The result is a long-format
tibble linking every `(representative_id, member_id)` pair to its source
`file_name`, making it straightforward to compare cluster membership
across multiple proteomes or organisms. Optionally, per-pair alignment
statistics from
[`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
can be incorporated without re-running `diamond deepclust`.

## Usage

``` r
deepclust_annotate(
  input_file,
  seq_type = "protein",
  format = "fasta",
  delete_corrupt_cds = TRUE,
  path = NULL,
  comp_cores = 1,
  mutual_cover = 80,
  approx_id = NULL,
  eval = NULL,
  deepclust_params = NULL,
  quiet = TRUE,
  cluster_table = NULL,
  realign = NULL
)
```

## Arguments

- input_file:

  a character string, a character vector of file paths, or a *named
  list* mapping custom labels to file paths. When a named list is
  supplied (e.g.
  `list("thal" = "/path/to/thal.fasta", "lyra" = "/path/to/lyra.fasta")`),
  the list *names* are used as the `file_name` column in the output
  instead of the base names of the files. This is useful when the file
  names are long or uninformative and you want cleaner labels in
  downstream analyses. Files are passed directly to
  [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  for clustering, and are also scanned individually to build a
  sequence-ID-to-file mapping (the "sequence library"). Ignored when
  `cluster_table` is supplied (the sequence library is still built from
  these paths in that case).

- seq_type:

  a character string specifying the sequence type stored in the input
  file(s). Options are: `"cds"`, `"protein"`, or `"dna"`. Default is
  `seq_type = "protein"`.

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

  a numeric value (percentage) specifying the minimum mutual coverage
  between a cluster member and its representative. Passed to
  `--mutual-cover`. Default is `mutual_cover = 80`.

- approx_id:

  a numeric value (percentage) specifying the minimum approximate
  identity for clustering. Passed to `--approx-id`. Default is `NULL`
  (DIAMOND2 default applies).

- eval:

  a numeric value specifying the maximum E-value. Passed to `--evalue`.
  Default is `NULL` (DIAMOND2 default of 0.001 applies).

- deepclust_params:

  a character string of additional DIAMOND2 deepclust parameters to
  append to the `deepclust` call. Default is `NULL`.

- quiet:

  a logical value indicating whether DIAMOND2 should run in quiet mode.
  Default is `quiet = TRUE`.

- cluster_table:

  an optional `tibble` (or `data.frame`) with columns
  `representative_id` and `member_id`, as returned by a previous call to
  [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
  or `deepclust_annotate`. When supplied, `diamond deepclust` is *not*
  re-run; the provided table is used directly as the clustering result.
  This is particularly useful when combined with `realign` to annotate
  an existing clustering with alignment statistics without incurring the
  cost of a second deepclust run. Default is `cluster_table = NULL`.

- realign:

  an optional `tibble` returned by
  [`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md).
  When supplied, the alignment statistics (`approx_pident`, `evalue`,
  `bitscore`, `qcovhsp`, `scovhsp`) are left-joined onto the cluster
  profile by `(representative_id, member_id)`, enriching each row with
  per-pair alignment information. Combine with `cluster_table` to
  annotate a previously computed clustering without re-running
  `diamond deepclust`. Default is `realign = NULL`.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
three columns (plus five additional alignment columns when `realign` is
supplied):

- `representative_id` — accession of the cluster representative.

- `member_id` — accession of the cluster member.

- `file_name` — base name of the source file from which the member
  sequence originated.

- `approx_pident` — (only with `realign`) approximate percentage of
  identical matches between representative and member.

- `evalue` — (only with `realign`) expect value.

- `bitscore` — (only with `realign`) bit score.

- `qcovhsp` — (only with `realign`) query (representative) coverage per
  HSP.

- `scovhsp` — (only with `realign`) subject (member) coverage per HSP.

## Details

The function performs up to four steps:

1.  **Build sequence library.** Each input file is read with
    [`seqinr::read.fasta`](https://rdrr.io/pkg/seqinr/man/read.fasta.html)
    to extract the sequence IDs (FASTA header names, i.e. the
    `member_id` values that DIAMOND2 will produce). These are collected
    into a two-column tibble of `member_id` and `file_name` (the base
    name of the source file).

2.  **Run deepclust** (skipped when `cluster_table` is supplied).
    [`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md)
    is called with all supplied `input_file` paths and the specified
    clustering parameters.

3.  **Join.** The sequence library is left-joined onto the deepclust
    output on `member_id`, annotating every cluster member with its
    source file.

4.  **Join realign** (only when `realign` is supplied). The alignment
    statistics from
    [`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md)
    are left-joined onto the profile by
    `(representative_id, member_id)`.

The resulting long-format tibble can be pivoted to a wide
presence/absence matrix (one row per cluster, one column per
organism/file) using
[`tidyr::pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html).

## References

Buchfink, B., Reuter, K., & Drost, H. G. (2021) "Sensitive protein
alignments at tree-of-life scale using DIAMOND." Nature methods, 18(4),
366-368.

https://github.com/bbuchfink/diamond/wiki

## See also

[`deepclust`](https://drostlab.github.io/orthologr/reference/deepclust.md),
[`deepclust_realign`](https://drostlab.github.io/orthologr/reference/deepclust_realign.md),
[`diamond`](https://drostlab.github.io/orthologr/reference/diamond.md),
[`set_diamond`](https://drostlab.github.io/orthologr/reference/set_diamond.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
if (FALSE) { # \dontrun{
# Profile two proteomes with default settings
profile <- deepclust_annotate(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  )
)

# Use a named list to control the file_name column labels
profile <- deepclust_annotate(
  input_file = list(
    "thal" = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    "lyra" = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  )
)
profile

# Pivot to a presence/absence matrix
library(tidyr)
library(dplyr)
profile |>
  distinct(representative_id, file_name) |>
  mutate(present = 1L) |>
  pivot_wider(
    names_from  = file_name,
    values_from = present,
    values_fill = 0L
  )

# Stricter clustering with custom thresholds
deepclust_annotate(
  input_file = c(
    system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
    system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
  ),
  approx_id    = 50,
  mutual_cover = 90,
  comp_cores   = 4
)

# Enrich an existing cluster profile with alignment stats from
# deepclust_realign, without re-running diamond deepclust
input_files <- c(
  system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
  system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr')
)

clusters      <- deepclust(input_file = input_files)
realign_stats <- deepclust_realign(input_file = input_files,
                                   clusters   = clusters)

profile_enriched <- deepclust_annotate(
  input_file    = input_files,
  cluster_table = clusters,
  realign       = realign_stats
)
} # }
```
