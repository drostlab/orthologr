# Sort dNdS Values Into Divergence Strata

This function takes a data.table returned by dNdS and sorts the
corresponding dNdS value into divergence strata (10-quantile or deciles
by default).

## Usage

``` r
divergence_map(dNdS_tbl, subject.id = FALSE, n_quantile = 10)
```

## Arguments

- dNdS_tbl:

  a data.table object returned by
  [`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md).

- subject.id:

  a logical value indicating whether `query_id` AND `subject_id` should
  be returned.

- n_quantile:

  a numeric value specifying the number of quantiles that should be
  returned.

## Value

a data.table storing a standard divergence map.

## Details

Divergence Strata are typically decile (10-quantile) values of
corresponding
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md) values.
The [`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)
function returns dNdS values for orthologous genes of a query species
(versus subject species). These dNdS values are then sorted into deciles
and each orthologous protein coding gene of the query species receives a
corresponding decile value instead of the initial dNdS value.

This allows a better comparison between Phylostrata and Divergence
Strata (for more details see package: myTAI).

It is also possible to specify other number of quantile (N-quantile)
value such as quintile (5-quantile) rather than decile values.

## See also

[`divergence_stratigraphy`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{

# get a divergence map of example sequences
dNdS_tbl <- dNdS( 
              query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
              subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
              ortho_detection = "RBH", 
              aa_aln_type     = "pairwise",
              aa_aln_tool     = "NW", 
              codon_aln_tool  = "pal2nal", 
              dnds_est.method = "Comeron", 
              comp_cores      = 1 )
                 
divMap <- divergence_map( dNdS_tbl )                       

# in case you want the subject_id as well, you can set
# the argument subject.id = TRUE
divMap <- divergence_map( dNdS_tbl = dNdS_tbl, subject.id = TRUE)   
            
 # in case you want a divergence map with divergence stratum as quintile (5-quantile) values 
 # or any other N-quantile values rather than the default decile (10-quantile) values,
 # you can specify this with n_quantile.
divMap <- divergence_map( dNdS_tbl = dNdS_tbl, subject.id = TRUE, n_quantile = 5) 

} # }
```
