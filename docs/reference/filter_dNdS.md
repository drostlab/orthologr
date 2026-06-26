# Filter dNdS values

This function takes the output data.table returned by
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md) and
filters the output by the following criteria:

1\) all dN values having an NA value are omitted

2\) all dS values having an NA value are omitted

3\) all dNdS values \>= the specified `dnds.threshold` are omitted

## Usage

``` r
filter_dNdS(dNdS_tbl, dnds.threshold = 2)
```

## Arguments

- dNdS_tbl:

  a `data.table` returned by
  [`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md).

- dnds.threshold:

  a numeric value specifying the dnds threshold for genes that shall be
  retained. Hence all genes having a dNdS value \<= `dnds.threshold` are
  retained. Default is `dnds.threshold` = 2.

## Details

The dNdS ratio quantifies the selection pressure acting on a given
protein sequence.

It is proposed that:

1\) dNdS values \< 1 reflect stabilizing selection of a protein

2\) dNdS value = 1 reflect neutral selection of a protein

3\) dNdS values \> 1 reflect variational selection of a protein

Now an assumption must be generated to allow for dNdS filtering. Given a
`dnds.threshold` all value above this threshold are omitted from the
dataset.

## See also

[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md),
[`divergence_stratigraphy`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{

filter_dNdS( dNdS( query_file = 
                   system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                   subject_file = 
                   system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
                   ortho_detection = "RBH", 
                   aa_aln_type = "pairwise",
                   aa_aln_tool = "NW",
                   codon_aln_tool = "pal2nal", 
                   dnds_est.method = "Comeron", 
                   comp_cores = 1), 
             dnds.threshold = 2)

# a small example using clustalw
filter_dNdS( dNdS( query_file      = 
                   system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
                   subject_file    = 
                   system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
                   ortho_detection = "RBH", 
                   aa_aln_type     = "pairwise",
                   aa_aln_tool     = "NW",
                   codon_aln_tool  = "pal2nal", 
                   dnds_est.method = "Comeron", 
                   comp_cores      = 1) , 
               dnds.threshold  = 2)
                  
} # }
```
