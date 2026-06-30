# Retrieve a core set of orthologous lncRNAs from the pairwise lncRNA orthologs map

Given a lnc_map table generated with
[`map_generator_lnc`](https://drostlab.github.io/orthologr/reference/map_generator_lnc.md),
this function will determine a core set of lncRNA orthologs that are
shared between all species.

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
  phylogenetic/taxonomic distance from the query species. The species
  names must match with the species names present in the `lnc_map`.

## Author

Hajk-Georg Drost
