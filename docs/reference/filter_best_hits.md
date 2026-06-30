# Helper function to select best BLAST hit based on minimum evalue

Helper function to select best BLAST hit based on minimum evalue

## Usage

``` r
filter_best_hits(x)
```

## Arguments

- x:

  a tibble storing gene locus ids, splice variant ids, and blast output
  for filtering

## Details

Extracts the best hit based on lowest e-value and longest splice variant
when e-values are equal

## Author

Hajk-Georg Drost
