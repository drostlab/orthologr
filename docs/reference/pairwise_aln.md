# Compute Pairwise Alignments

This function takes a FASTA file containing two DNA or amino acid
sequences that shall be aligned and computes a paiwise alignment using a
defined alignment method.

## Usage

``` r
pairwise_aln(
  file,
  tool = "NW",
  seq_type = "protein",
  get_aln = FALSE,
  pairwise_aln_name = NULL,
  path = NULL,
  quiet = FALSE,
  store_locally = FALSE
)
```

## Arguments

- file:

  a character string specifying the path to the file storing the
  sequences in FASTA format.

- tool:

  a character string specifying the program/algorithm that should be
  used: `tool` = `"NW"`.

- seq_type:

  a character string specifying the sequence type stored within the
  given FASTA file. Options are `seq_type` = `"protein"`, `seq_type` =
  `"cds"`, `seq_type` = `"dna"`. Default is `seq_type` = `"protein"`.

- get_aln:

  a logical value indicating whether the produced alignment should be
  returned.

- pairwise_aln_name:

  a character string specifying the name of the stored alignment file.
  Default is `pairwise_aln_name` = `NULL` denoting a default name:
  'toolname_seq_type.aln' .

- path:

  a character string specifying the path to the pairwise alignment
  program (in case you don't use the default path).

- quiet:

  a logical value specifying whether a successful interface call shall
  be printed out.

- store_locally:

  a logical value indicating whether or not alignment files shall be
  stored locally rather than in
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

## Value

In case the argument `get_aln` is set `TRUE`, an object of class
alignment of the seqinr package is returned.

## Details

This function provides an interface between R and common pairwise
alignment computation methods.

The current version of this function computes pairwise alignments based
on the
[`pairwiseAlignment`](https://rdrr.io/pkg/pwalign/man/pairwiseAlignment.html)
function implemented in the pwalign package.

The default pairwise alignment method is based on the *Needleman-Wunsch
Algorithm*.

## Author

Sarah Scharfenberg and Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{        
                                            
# Needleman-Wunsch Example:  

# in case Biostrings works properly
pairwise_aln( file     = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
              tool     = "NW", 
              get_aln  = TRUE, 
              seq_type = "protein")

                                                   
                                                                                                                                                         
} # }
```
