# Estimate the DNA distance between promotor sequences

Given two `fasta` files storing the promotor sequence of protein coding
genes, this function estimates the pairwise DNA distance between these
promotor sequences.

## Usage

``` r
promotor_divergence_estimation(query, subject, model = "K80")
```

## Arguments

- query:

  file path to a query `fasta` file storing promotor sequences of
  interest (e.g. generated with
  [`extract_upstream_promotor_seqs`](https://rdrr.io/pkg/metablastr/man/extract_upstream_promotor_seqs.html)).

- subject:

  file path to a subject `fasta` file storing promotor sequences of
  interest (e.g. generated with
  [`extract_upstream_promotor_seqs`](https://rdrr.io/pkg/metablastr/man/extract_upstream_promotor_seqs.html)).

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

## Author

Hajk-Georg Drost
