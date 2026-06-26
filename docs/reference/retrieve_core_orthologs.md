# Retrieve a core set of orthologous gene loci from several pairwise ortholog tables

Given a ortho table generated with
[`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md),
this function will determine a core set of orthologs that are shared
between all species.

## Usage

``` r
retrieve_core_orthologs(ortho_tables, species_order)
```

## Arguments

- ortho_tables:

  a `ortho tables` that was generated with
  [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md).

- species_order:

  a character string specifying species names listed in the order of
  phylogenetic/taxonomic distance from the query species. The species
  names must match with the species names present in the `ortho_tables`.

## Author

Hajk-Georg Drost
