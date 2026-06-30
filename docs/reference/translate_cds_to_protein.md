# Translate coding sequences into amino acid sequences

A helper function that takes a `fasta` file storing coding sequences as
input and translates these coding sequences into amino acid sequences
storing them as `fasta` output file.

## Usage

``` r
translate_cds_to_protein(input_file, output_file, delete_corrupt_cds = FALSE)
```

## Arguments

- input_file:

  file path to `fasta` file storing coding sequences (DNA).

- output_file:

  name or file path in which translated amino acid sequences (AA) shall
  be stored as `fasta` file.

- delete_corrupt_cds:

  delete_corrupt_cds a logical value indicating whether sequences with
  corrupt base triplets should be removed from the input `file`. This is
  the case when the length of coding sequences cannot be divided by 3
  and thus the coding sequence contains at least one corrupt base
  triplet.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# install.packages("biomartr")
# download coding sequences of Arabidopsis thaliana
Ath_file_path <- biomartr::getCDS(db = "refseq",
                     organism = "Arabidopsis thaliana",
                     gunzip = TRUE)
# translate coding sequences into amino acid sequences
translate_cds_to_protein(Ath_file_path, "Ath_aa_seqs.fasta")
} # }
```
