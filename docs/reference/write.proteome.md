# Save a proteome in fasta format

This function takes a proteome (in data.table notation) as returned by
[`read.proteome`](https://drostlab.github.io/orthologr/reference/read.proteome.md)
and saves it to your hard drive.

## Usage

``` r
write.proteome(proteome, file.name, nbchar = 80, open = "w", as.string = TRUE)
```

## Arguments

- proteome:

  a data.table object storing amino acid sequences in the first column
  and gene ids in the second column. See
  [`read.proteome`](https://drostlab.github.io/orthologr/reference/read.proteome.md)
  for details.

- file.name:

  a character string specifying the name of the fasta output file.

- nbchar:

  number of characters per line.

- open:

  mode to open the output file, you can choose from "w" to write into a
  new file, or "a" to append at the end of an already existing file.

- as.string:

  when set to TRUE sequences are in the form of strings instead of
  vectors of single characters.

## Author

Hajk-Georg Drost

## Examples

``` r

if (FALSE) { # \dontrun{

# read an example proteome
Ath.proteome <- read.proteome(
                      system.file('seqs/ortho_thal_aa.fasta', package = 'orthologr'), 
                      format = "fasta")

# use the write.proteome function to store it on your hard drive
write.proteome(Ath.proteome, "test_proteome.fasta")

} # }
```
