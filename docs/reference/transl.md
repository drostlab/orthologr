# Translate DNA to Amino Acids

This function takes CDS sequence as string as input an returns the
corresponding amino acid sequence as string.

## Usage

``` r
transl(sequence)
```

## Arguments

- sequence:

  a character string specifying CDS sequence of interest.

## Value

A character string specifying the corresponding amino acid sequence.

## References

[`translate`](https://rdrr.io/pkg/seqinr/man/translate.html)

## Author

Hajk-Georg Drost and Sarah Scharfenberg

## Examples

``` r

# an example DNA sequence
DNA <- c("ACCGGTTTAAAGGCGTTA")

# translating DNA to a protein sequence
transl(DNA)
#> [1] "TGLKAL"
```
