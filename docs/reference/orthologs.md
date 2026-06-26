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
  ortho_detection = "RBH",
  delete_corrupt_cds = FALSE,
  cdd.path = NULL,
  path = NULL,
  add_params = NULL,
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
  interest (subject organisms). Different orthology inference methods
  can detect orthologs using multiple subject organisms, e.g.
  "OrthoMCL", and "PO" (ProteinOrtho).

- seq_type:

  a character string specifying the sequence type stored in the input
  file. Options are are: "cds", "protein", or "dna". In case of "cds",
  sequence are translated to protein sequences, in case of "dna", cds
  prediction is performed on the corresponding sequences which
  subsequently are translated to protein sequences. Default is
  `seq_type` = "protein".

- outgroup_file:

  a character string specifying the paths to the sequence files of
  interest (outgroup organisms). This argument is only used by
  `InParanoid`.

- eval:

  a numeric value specifying the E-Value cutoff for BLAST hit detection.

- format:

  a character string specifying the file format of the sequence file,
  e.g. "fasta", "gbk". Default is "fasta".

- ortho_detection:

  a character string specifying the orthology inference method that
  shall be performed to detect orthologous genes. Default is
  `ortho_detection` = "RBH" (BLAST reciprocal best hit). Further methods
  are: "RBH" (BLAST reciprocal best hit), "PO" (ProteinOrtho), and
  "OrthoMCL.

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- cdd.path:

  path to the cdd database folder (specify when using `ortho_detection`
  = `"DELTA"`).

- path:

  a character string specifying the path to the corresponding orthology
  inference tool. For "BH" and "RBH": path to BLAST, "PO": path to
  ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL.

- add_params:

  a character string specifying additional parameters that shall be
  handed to the orthology inference method (tool). Default is
  `add_params` = `NULL`.

- comp_cores:

  a numeric value specifying the number of cores to be used for
  multicore computations.

- quiet:

  a logical value specifying whether a successful interface call shall
  be printed out.

- clean_folders:

  a boolean value spefiying whether all internall folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

## Value

A data.table storing the query_ids of orthologous genes in the first
column, the subject_ids of orthologous genes in the second column and
the amino acid sequences in the third column.

## Details

This function takes sequence files of a query organism and a subject
organism and performs orthology inference using a defined orthology
inference method to dectect orthologous genes.

The following interfaces are implemented in the `orthologs` function:

BLAST based methods:

- BLAST best hit (BH)

- BLAST reciprocal best hit (RBH)

- DELTA-BLAST reciprocal best hit (DELTA)

## See also

[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{

### BLAST Best Hit

# perform orthology inference using BLAST best hit
# and fasta sequence files storing protein sequences
orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein", 
          ortho_detection = "BH")


### BLAST Reciprocal Best Hit

# perform orthology inference using BLAST reciprocal best hit
# and fasta sequence files storing protein sequences
orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein", 
          ortho_detection = "RBH")
          
          
# multicore version          
orthologs(query_file      = system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
          subject_files   = system.file('seqs/ortho_lyra_aa.fasta', package = 'orthologr'),
          seq_type        = "protein", 
          ortho_detection = "RBH", 
          comp_cores      = 2)          
} # }
```
