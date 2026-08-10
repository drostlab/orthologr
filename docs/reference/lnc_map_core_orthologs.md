# Retrieve a core set of orthologous lncRNAs (deprecated)

This function is deprecated and will be removed in a future release.
Please use
[`retrieve_core_orthologs`](https://drostlab.github.io/orthologr/reference/retrieve_core_orthologs.md)
instead, which now handles both protein-coding gene tables and lncRNA
maps automatically.

## Usage

``` r
lnc_map_core_orthologs(lnc_map, species_order)
```

## Arguments

- lnc_map:

  a `lnc_map` that was generated with
  [`map_generator_lnc`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md).

- species_order:

  a character string specifying species names listed in the order of
  phylogenetic/taxonomic distance from the query species.

## Author

Hajk-Georg Drost
