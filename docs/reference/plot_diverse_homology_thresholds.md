# Diverse line plots visualizing the number of pairwise orthologs within a ortho table generated with [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md) based on different sets of homology thresholds.

Given a ortho table generated with
[`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md),
this function will visualize the number of pairwise orthologs inferred
between a reference species A vs a set of subject species B_1, B_2,
...,B_N using a diverse range of homology thresholds to enable an
analytical decision for chosing homology thresholds to define
orthologous genes.

## Usage

``` r
plot_diverse_homology_thresholds(
  ortho_tables,
  species_order,
  xlab = "Subject Species",
  ylab = "Number of reciprocal best hit orthologs",
  title = ""
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

## Author

Hajk-Georg Drost
