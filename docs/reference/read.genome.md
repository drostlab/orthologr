# Read the genome of a given organism

This function reads an organism specific genome stored in a defined file
format.

## Usage

``` r
read.genome(file, format, ...)
```

## Arguments

- file:

  a character string specifying the path to the file storing the genome.

- format:

  a character string specifying the file format used to store the
  genome, e.g. "fasta", "fastq".

- ...:

  additional arguments that are used by the
  [`readDNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
  function.

## Value

A data.table storing the gene id in the first column and the
corresponding sequence as string in the second column.

## Details

The `read.genome` function takes a string specifying the path to the
genome file of interest as first argument.

It is possible to read in different genome file standards such as
*fasta* or *fastq*. Genomes stored in fasta files can be downloaded from
https://ensemblgenomes.org, https://www.ncbi.nlm.nih.gov/genome, or
https://www.ebi.ac.uk, etc.

## Author

Hajk-Georg Drost and Sarah Scharfenberg

## Examples

``` r
if (FALSE) { # \dontrun{
# reading a genome stored in a fasta file
Ath.genome <- read.genome(system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                           format = "fasta")
} # }
```
