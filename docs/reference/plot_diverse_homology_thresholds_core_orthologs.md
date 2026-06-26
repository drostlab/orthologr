# Diverse line plots visualizing the number of core orthologs within a ortho table generated with [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md) based on different sets of homology thresholds.

Given a ortho table generated with
[`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md),
this function will visualize the number of core orthologs using a
diverse range of homology thresholds to enable an analytical decision
for chosing homology thresholds to define orthologous genes.

## Usage

``` r
plot_diverse_homology_thresholds_core_orthologs(
  ortho_tables,
  species_order,
  type = "both",
  xlab = "Parameter set for selecting core orthologs",
  ylab = "Number of corresponding core orthologs",
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

- type:

  type of core orthologs that shall be visualized:

  - `type = "gene_locus"`: visualizes the number of core orthologs for a
    gene locus (hence, different splice variants can represent the same
    orthologous gene locus relationship)

  - `type = "both"`: visualizes the difference in numbers of core
    orthologs for a gene locus versus orthologous splice variant (thus,
    only the same splice variant can represent the same orthologous gene
    locus relationship)

- xlab:

  x-axis label.

- ylab:

  y-axis label.

- title:

  plot title.

## Author

Hajk-Georg Drost
