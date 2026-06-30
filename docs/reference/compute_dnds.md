# Compute dNdS Values For A Given Pairwise Alignment

This function takes a vector containing the query amino acid sequence,
subject amino acid sequence, query CDS sequence, and subject CDS
sequence and then runs the following pipieline:

- 1\) Multiple-Alignment of query amino acid sequence and subject amino
  acid sequence

- 2\) Codon-Alignment of the amino acid alignment returned by 1) and
  query CDS sequence + subject CDS sequence

- 3\) dNdS estimation of the codon alignment returned by 2)

## Usage

``` r
compute_dnds(
  complete_tbl,
  aa_aln_type = "multiple",
  aa_aln_tool = "clustalw",
  aa_aln_path = NULL,
  aa_aln_params = NULL,
  codon_aln_tool = "pal2nal",
  dnds_est.method = "YN",
  store_locally = FALSE,
  kaks_calc_path = NULL,
  quiet = FALSE,
  comp_cores = 1,
  clean_folders = FALSE
)
```

## Arguments

- complete_tbl:

  a data.table object storing the query_id, subject_id, query_cds
  (sequence), subject_cds (sequence), query_aa (sequence), and
  subject_aa (sequence) of the organisms that shall be compared.

- aa_aln_type:

  a character string specifying the amino acid alignement type:
  `aa_aln_type` = "multiple" or `aa_aln_type` = `"pairwise"`. Default is
  `aa_aln_type` = "multiple".

- aa_aln_tool:

  a character string specifying the multiple alignment tool that shall
  be used for pairwise protein alignments.

- aa_aln_path:

  a character string specifying the path to the corresponding multiple
  alignment tool.

- aa_aln_params:

  a character string specifying additional parameters that shall be
  passed to the multiple alignment system call.

- codon_aln_tool:

  a character string specifying the codon alignment tool that shall be
  used for codon alignments. Default is `codon_aln_tool` = `"pal2nal"`.

- dnds_est.method:

  a character string specifying the dNdS estimation method, e.g.
  "Comeron","Li", "YN", etc. See Details for all options.

- store_locally:

  a logical value indicating whether or not alignment files shall be
  stored locally rather than in
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- kaks_calc_path:

  a character string specifying the execution path to KaKs_Calculator.
  Default is `kaks_calc_path` = `NULL` (meaning that KaKs_Calculator is
  stored and executable in your default `PATH`).

- quiet:

  a logical value specifying whether a successful interface call shall
  be printed out.

- comp_cores:

  a numeric value specifying the number of cores that shall be used to
  perform parallel computations on a multicore machine.

- clean_folders:

  a boolean value spefiying whether all internall folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

## Details

This function takes the amino acid and CDS sequences two orthologous
genes and writes the corresponding amino acid and CDS sequences as fasta
file into the internal folder environment. The resulting fasta files
(two files) store the amino acid sequence of the query_id and subject_id
(file one) and the CDS sequence of the query_id and subject_id (file
two). These fasta files are then used to pass through the following
pipeline:

1\) Multiple-Alignment or Pairwise-Alignment of query amino acid
sequence and subject amino acid sequence

2\) Codon-Alignment of the amino acid alignment returned by 1) and query
CDS sequence + subject CDS sequence

3\) dNdS estimation of the codon alignment returned by 2)

## References

<https://www.r-bloggers.com/2013/08/the-wonders-of-foreach/>

## See also

[`multi_aln`](https://drostlab.github.io/orthologr/reference/multi_aln.md),
[`substitutionrate`](https://drostlab.github.io/orthologr/reference/substitutionrate.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)

## Author

Hajk-Georg Drost and Sarah Scharfenberg
