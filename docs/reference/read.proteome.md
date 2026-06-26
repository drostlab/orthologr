# Read the proteome of a given organism

This function reads an organism specific proteome stored in a defined
file format.

## Usage

``` r
read.proteome(file, format, ...)
```

## Arguments

- file:

  a character string specifying the path to the file storing the
  proteome.

- format:

  a character string specifying the file format used to store the
  proteome, e.g. "fasta", "fastq".

- ...:

  additional arguments that are used by the
  [`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
  function.

## Value

A data.table storing the gene id in the first column and the
corresponding sequence as string in the second column.

## Details

The `read.proteome` function takes a string specifying the path to the
proteome file of interest as first argument.

It is possible to read in different proteome file standards such as
*fasta* or *fastq*.

Proteomes stored in fasta files can be downloaded from
https://www.ebi.ac.uk/reference_proteomes, etc.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# reading a proteome stored in a fasta file
Ath.proteome <- read.proteome(system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'),
                               format = "fasta")
} # }
```
