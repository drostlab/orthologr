# Helper function to extract a core set of orthologs

Helper function to extract a core set of orthologs shared across all
subject species. Handles both protein-coding gene tables (species column
`subject_species`) and lncRNA maps (species column `species`). The
appropriate column is detected automatically from the input.

## Usage

``` r
filter_core_set(x, order_species)
```

## Arguments

- x:

  input data in `data.frame` or `tibble` format.

- order_species:

  a character vector containing the scientific names of the organisms of
  interest ordered according to their phylogenetic distance to their
  reference species.

## Examples

``` r
if (FALSE) { # \dontrun{
# Protein-coding table: all target species present -> returns x unchanged
x_complete <- tibble::tibble(
  subject_species     = c("Arabidopsis_lyrata", "Brassica_rapa"),
  query_id            = c("AT1G01010.1", "AT1G01010.1"),
  query_gene_locus_id = c("AT1G01010", "AT1G01010"),
  dN = c(0.05, 0.08), dS = c(0.12, 0.15), dNdS = c(0.42, 0.53)
)
filter_core_set(x_complete, order_species = c("Arabidopsis_lyrata", "Brassica_rapa"))

# Protein-coding table: one species missing -> returns a schema-preserving NA row
x_incomplete <- tibble::tibble(
  subject_species     = "Arabidopsis_lyrata",
  query_id            = "AT1G01010.1",
  query_gene_locus_id = "AT1G01010",
  dN = 0.05, dS = 0.12, dNdS = 0.42
)
filter_core_set(x_incomplete, order_species = c("Arabidopsis_lyrata", "Brassica_rapa"))
} # }
```
