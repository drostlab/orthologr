# Read the CDS of a given organism

This function reads an organism specific CDS stored in a defined file
format.

## Usage

``` r
read.cds(file, format, delete_corrupt_cds = TRUE, ...)
```

## Arguments

- file:

  a character string specifying the path to the file storing the CDS.

- format:

  a character string specifying the file format used to store the CDS,
  e.g. "fasta", "fatsq".

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.

- ...:

  additional arguments that are used by the
  [`readDNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
  function.

## Value

A data.table storing the gene id in the first column and the
corresponding sequence as string in the second column.

## Details

The `read.cds` function takes a string specifying the path to the cds
file of interest as first argument.

It is possible to read in different proteome file standards such as
*fasta* or *fastq*.

CDS stored in fasta files can be downloaded from
https://www.ensembl.org/info/data/ftp/index.html, etc.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# reading a cds file stored in fasta format
Ath.cds <- read.cds(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                    format = "fasta")
} # }
```
