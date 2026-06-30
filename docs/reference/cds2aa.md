# Translate CDS file to Amino Acids file

This function takes a `data.table` object as input and translates the
cds sequences stored as `seqs` column into the corresponding amino acid
sequence.

## Usage

``` r
cds2aa(dt, delete_corrupt_cds = TRUE)
```

## Arguments

- dt:

  a `data.table` object storing cds sequences in a column named `seqs`.

- delete_corrupt_cds:

  a logical value indicating whether sequences with corrupt base
  triplets should be removed from the input `file`. This is the case
  when the length of coding sequences cannot be divided by 3 and thus
  the coding sequence contains at least one corrupt base triplet.
