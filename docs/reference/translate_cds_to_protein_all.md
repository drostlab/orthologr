# Translate coding sequences into amino acid sequences for multiple files

A helper function that takes `fasta` files storing coding sequences as
input and translates these coding sequences into amino acid sequences
storing them as `fasta` output files.

## Usage

``` r
translate_cds_to_protein_all(
  input_folder,
  output_folder,
  delete_corrupt_cds = FALSE
)
```

## Arguments

- input_folder:

  file path to folder storing the coding sequences `fasta` files.

- output_folder:

  name or file path to a folder that that shall be generated to store
  the output `fasta` files.

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
# download coding sequences of Arabidopsis thaliana, Arabidopsis lyrata, and Capsella rubella
org_list <- c("Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella")
biomartr::getCDSSet(db = "refseq",
             organism = org_list,
             gunzip = TRUE,
             path = "cds_examples")
# translate coding sequences into amino acid sequences
translate_cds_to_protein_all(input_folder = "cds_examples", 
                             output_folder = "translated_seqs",
                             delete_corrupt_cds = FALSE)
} # }
```
