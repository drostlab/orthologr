# Main Orthology Inference Function

This function takes nucleotide or protein sequences for a set of
organisms and performs orthology inference to detect orthologous genes
within the given organisms based on selected orthology inference
programs.

## Usage

``` r
orthologs(
  query_file,
  subject_files,
  seq_type = "protein",
  outgroup_file = NULL,
  eval = "1E-5",
  format = "fasta",
  ortho_detection = "DIAMOND_RBH",
  sensitivity_mode = "fast",
  max.target.seqs = 10000,
  delete_corrupt_cds = FALSE,
  cdd.path = NULL,
  path = NULL,
  add_params = NULL,
  diamond_params = NULL,
  comp_cores = 1,
  quiet = FALSE,
  clean_folders = FALSE
)
```

## Arguments

- query_file:

  a character string specifying the path to the sequence file of
  interest (query organism).

- subject_files:

  a character string specifying the paths to the sequence files of
  interest (subject organisms).

- seq_type:

  a character string specifying the sequence type stored in the input
  file. Options are: "cds", "protein", or "dna". In case of "cds",
  sequences are translated to protein sequences, in case of "dna", CDS
  prediction is performed on the corresponding sequences which are
  subsequently translated to protein sequences. Default is `seq_type` =
  "protein".

- outgroup_file:

  currently unused; reserved for future use.

- eval:

  a numeric value specifying the E-Value cutoff for hit detection.

- format:

  a character string specifying the file format of the sequence file,
  e.g. "fasta". Default is "fasta".

- ortho_detection:

  a character string specifying the orthology inference method to use.
  Default is `ortho_detection` = `"DIAMOND_RBH"` (DIAMOND2 reciprocal
  best hit). Available methods:

  - `"DIAMOND_RBH"`: DIAMOND2 reciprocal best hit (default; fast,
    recommended)

  - `"DIAMOND_BH"`: DIAMOND2 best hit

  - `"RBH"`: BLAST reciprocal best hit

  - `"BH"`: BLAST best hit

  - `"Orthofinder2"`: OrthoFinder2 (currently under development)

- sensitivity_mode:

  a character string specifying the DIAMOND2 sensitivity level. Only
  used when `ortho_detection` is `"DIAMOND_RBH"` or `"DIAMOND_BH"`.
  Options: `"faster"`, `"fast"` (default), `"mid-sensitive"`,
  `"sensitive"`, `"more-sensitive"`, `"very-sensitive"`,
  `"ultra-sensitive"`.

- max.target.seqs:

  a numeric value specifying the number of aligned sequences to keep
  (DIAMOND2 only). Default is `10000`.

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. Default is `FALSE`.

- cdd.path:

  currently unused; reserved for future use.

- path:

  a character string specifying the path to the orthology inference tool
  executable.

- add_params:

  a character string specifying additional parameters passed to the
  orthology inference tool. For DIAMOND2 methods, use `diamond_params`
  instead. Default is `NULL`.

- diamond_params:

  a character string specifying additional parameters passed to
  DIAMOND2. Only used when `ortho_detection` is `"DIAMOND_RBH"` or
  `"DIAMOND_BH"`. Default is `NULL`.

- comp_cores:

  a numeric value specifying the number of cores to be used for
  multicore computations.

- quiet:

  a logical value specifying whether a successful interface call shall
  be printed out.

- clean_folders:

  a boolean value specifying whether all internal folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

## Value

A data.table storing the query_ids of orthologous genes in the first
column, the subject_ids of orthologous genes in the second column and
the amino acid sequences in the third column.

## Details

This function takes sequence files of a query organism and a subject
organism and performs orthology inference using a defined orthology
inference method to detect orthologous genes.

The following interfaces are implemented in the `orthologs` function:

DIAMOND2 based methods (recommended):

- DIAMOND2 reciprocal best hit (`"DIAMOND_RBH"`) — default

- DIAMOND2 best hit (`"DIAMOND_BH"`)

BLAST based methods:

- BLAST reciprocal best hit (`"RBH"`)

- BLAST best hit (`"BH"`)

## See also

[`diamond_rec`](https://drostlab.github.io/orthologr/reference/diamond_rec.md),
[`diamond_best`](https://drostlab.github.io/orthologr/reference/diamond_best.md),
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{

### DIAMOND2 Reciprocal Best Hit (default)

orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein")


### DIAMOND2 Best Hit

orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "DIAMOND_BH")


### BLAST Reciprocal Best Hit

orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "RBH")


### BLAST Best Hit

orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein",
          ortho_detection = "BH")
} # }
```
