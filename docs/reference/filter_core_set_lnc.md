# Helper function to extract a core set of orthologous lncRNAs (deprecated)

This function is deprecated and will be removed in a future release.
Please use
[`filter_core_set`](https://drostlab.github.io/orthologr/reference/filter_core_set.md)
instead, which now handles both protein-coding gene tables and lncRNA
maps automatically.

## Usage

``` r
filter_core_set_lnc(x, order_species)
```

## Arguments

- x:

  input data in `data.frame` or `tibble` format.

- order_species:

  a character vector containing the scientific names of the organisms of
  interest ordered according to their phylogenetic distance to their
  reference species.
