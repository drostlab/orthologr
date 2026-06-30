# Helper function for splice variant based `plot_diverse_homology_thresholds_core_orthologs`

Helper function for splice variant based
[`plot_diverse_homology_thresholds_core_orthologs`](https://drostlab.github.io/orthologr/reference/plot_diverse_homology_thresholds_core_orthologs.md).

## Usage

``` r
testCoreOrthoParamsSpliceVariant(
  ortho_tables,
  param,
  species_order,
  n_core_species
)
```

## Arguments

- ortho_tables:

  a `ortho tables` that was generated with
  [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md).

- param:

  parameters.

- species_order:

  a character string specifying species names listed in the order of
  phylogenetic/taxonomic distance from the query species.

- n_core_species:

  number of core species.
