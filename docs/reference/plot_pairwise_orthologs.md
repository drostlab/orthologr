# A line plot visualizing the number of pairwise orthologs within a ortho table generated with [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md)

Given a ortho table generated with
[`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md),
this function will visualize the number of pairwise orthologs inferred
between a reference species A vs a set of subject species B_1, B_2,
...,B_N.

## Usage

``` r
plot_pairwise_orthologs(
  ortho_tables,
  species_order,
  xlab = "Subject Species",
  ylab = "Number of reciprocal best hit orthologs",
  title = "",
  n_core_orthologs = NULL
)
```

## Arguments

- ortho_tables:

  a `ortho tables` that was generated with
  [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md).

- species_order:

  a character string specifying species names listed in the order of
  phylogenetic/taxonomic distance from the query species. The species
  names must match with the species names present in the `ortho_tables`.

- xlab:

  x-axis label.

- ylab:

  y-axis label.

- title:

  plot title.

- n_core_orthologs:

  number of core orthologs within the `ortho tables`. This number can be
  retrieved using
  [`retrieve_core_orthologs`](https://drostlab.github.io/orthologr/reference/retrieve_core_orthologs.md).

## Author

Hajk-Georg Drost
