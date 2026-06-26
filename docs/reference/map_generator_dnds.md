# Generate dNdS Maps Between a Query Organism and Multiple Subject Organisms

This function allows you to compute dNdS or DS Maps between a query
organism and a set subject organisms stored in the same folder. The
corresponding dNdS/DS Maps are then stored in an output folder.

## Usage

``` r
map_generator_dnds(
  query_file,
  subjects_folder,
  output_folder,
  eval = "1E-5",
  min_qry_coverage_hsp = 50,
  min_qry_perc_identity = 30,
  ortho_detection = "RBH",
  aa_aln_type = "pairwise",
  aa_aln_tool = "NW",
  codon_aln_tool = "pal2nal",
  dnds_est_method = "Comeron",
  comp_cores = 1,
  progress_bar = TRUE,
  sep = ";",
  ...
)
```

## Arguments

- query_file:

  a character string specifying the path to the CDS file of the query
  organism.

- subjects_folder:

  a character string specifying the path to the folder where CDS files
  of the subject organisms are stored.

- output_folder:

  a character string specifying the path to the folder where output
  dNdS/DS Maps should be stored.

- eval:

  a character string specifying the e-value for BLAST based Orthology
  Inference that is performed in the process of dNdS computations.
  Please use the scientific notation.

- min_qry_coverage_hsp:

  minimum `qcovhsp` (= query coverage of the HSP) of an orthologous hit
  (a value between 1 and 100).

- min_qry_perc_identity:

  minimum `perc_identity` (= percent sequence identity between query and
  selected HSP) of an orthologous hit (a value between 1 and 100).

- ortho_detection:

  a character string specifying the Orthology Inference method that
  shall be used to perform dNdS computations. Possible options are:
  `ortho_detection` = `"BH"` (BLAST best hit), `ortho_detection` =
  `"RBH"` (BLAST best reciprocal hit), etc.

- aa_aln_type:

  a character string specifying the amino acid alignement type:
  `aa_aln_type = "multiple"` or `aa_aln_type = "pairwise"`. Default is
  `aa_aln_type = "pairwise"`.

- aa_aln_tool:

  a character string specifying the program that should be used e.g.
  "clustalw".

- codon_aln_tool:

  a character string specifying the codon alignment tool that shall be
  used. Default is `codon_aln_tool = "pal2nal"`. Right now only
  "pal2nal" can be selected as codon alignment tool.

- dnds_est_method:

  a character string specifying the dNdS estimation method, e.g.
  "Comeron","Li", "YN", etc. See Details for all options.

- comp_cores:

  number of computing cores that shall be used to perform parallelized
  computations.

- progress_bar:

  should a progress bar be shown. Default is `progress_bar = TRUE`.

- sep:

  a file separator that is used to store maps as csv file.

- ...:

  additional parameters that shall be passed to
  [`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md).

## Details

Given a query organism and a set of subject organsisms that are stored
in the same folder, this function crawls through all subject organsism
and as a first step computes the pairwise dNdS Maps between query and
subject organsim and as a second step stores the corresponding Map in an
output folder.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
map_generator_dnds(
   query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
   subjects_folder = system.file('seqs/map_gen_example', package = 'orthologr'),
   aa_aln_type      = "pairwise",
   aa_aln_tool      = "NW", 
   codon_aln_tool   = "pal2nal", 
   dnds_est_method  = "Comeron",
   output_folder   = "orthologr_dnds_maps",
   quiet           = TRUE,
   comp_cores      = 1
)

} # }
```
