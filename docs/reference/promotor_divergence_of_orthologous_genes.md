# Compute promotor sequence divergence of orthologous genes

This function computes the promotor sequence divergences of orthologous
genes from a set of pairwise species comparisons. It allows to add this
promotor divergence information to an pre-computed `dNdS table`
generated with
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md) or
[`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md).

## Usage

``` r
promotor_divergence_of_orthologous_genes(
  promotor_folder,
  ortholog_tables_folder,
  model = "K80",
  ortholog_promotor_seq_output = NULL
)
```

## Arguments

- promotor_folder:

  a path to a folder containing promotor sequences in `fasta` format
  that were generated with
  [`extract_upstream_promotor_seqs`](https://rdrr.io/pkg/metablastr/man/extract_upstream_promotor_seqs.html).

- ortholog_tables_folder:

  a path to a folder containing dNdS tables generated with
  [`generate_ortholog_tables_all`](https://drostlab.github.io/orthologr/reference/generate_ortholog_tables_all.md)

- model:

  a model as specified in
  [`dist.dna`](https://rdrr.io/pkg/ape/man/dist.dna.html): a character
  string specifying the evolutionary model to be used - must be one of:

  - `K80` (the default)

  - `raw`

  - `N`

  - `TS`

  - `TV`

  - `JC69`

  - `F81`

  - `K81`

  - `F84`

  - `BH87`

  - `T92`

  - `TN93`

  - `GG95`

  - `logdet`

  - `paralin`

- ortholog_promotor_seq_output:

  a path or name to an output folder where orthologous promotor
  sequences shall be stored.

## Author

Hajk-Georg Drost
