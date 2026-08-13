# Infer orthologous lncRNAs between multiple species

Inference of orthologous lncRNAs between multiple species is performed
via pairwise BLAST (reciprocal) best hit comparisons. The corresponding
orthologous tables are then stored in an output folder.

## Usage

``` r
map_generator_lnc(
  query_file,
  subjects_folder,
  output_folder,
  task = "blastn",
  eval = "1E-5",
  ortho_detection = "RBH",
  max.target.seqs = 10000,
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  logical_connective = "AND",
  min_alig_length = NULL,
  comp_cores = 1,
  progress_bar = TRUE,
  sep = ";",
  path = NULL,
  ...
)
```

## Arguments

- query_file:

  a character string specifying the path to the lncRNAs file of the
  query organism in `fasta` format.

- subjects_folder:

  a character string specifying the path to the folder where lncRNAs
  files in `fasta` format of the subject organisms are stored.

- output_folder:

  a character string specifying the path to the folder where output
  orthologous tables should be stored.

- task:

  nucleotide search task option. Options are:

  - `task = "blastn"` : Standard nucleotide-nucleotide comparisons
    (default) - Traditional BLASTN requiring an exact match of 11.

  - `task = "blastn-short"` : Optimized nucleotide-nucleotide
    comparisons for query sequences shorter than 50 nucleotides.

  - `task = "dc-megablast"` : Discontiguous megablast used to find
    somewhat distant sequences.

  - `task = "megablast"` : Traditional megablast used to find very
    similar (e.g., intraspecies or closely related species) sequences.

  - `task = "rmblastn"`

- eval:

  a character string specifying the e-value for BLAST based orthology
  inference. Please use the scientific notation.

- ortho_detection:

  a character string specifying the Orthology Inference method that
  shall be used to perform dNdS computations. Possible options are:

  - `ortho_detection = "BH"`: BLAST best unidirectional hit

  - `ortho_detection = "RBH"`: BLAST best reciprocal hit

- max.target.seqs:

  a numeric value specifying the number of aligned sequences to keep.
  Please be aware that `max.target.seqs` selects best hits based on the
  database entry and not by the best e-value. See details here:
  https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166
  .

- min_qry_coverage_hsp:

  minimum `qcovhsp` (= query coverage of the HSP) of an orthologous hit
  (a value between 1 and 100).

- min_qry_perc_identity:

  minimum `perc_identity` (= percent sequence identity between query and
  selected HSP) of an orthologous hit (a value between 1 and 100).

- logical_connective:

  character representing logical connective (either "AND" or "OR") if
  `min_alig_length` is not `NULL` filtering is done on `min_alig_length`
  and/or `min_qry_perc_identity`

- min_alig_length:

  minimum `alig_length` (alignment length) to an orthologous hit (number
  of aligned nucleotides or amino acids depending on the input data)

- comp_cores:

  number of computing cores that shall be used to perform parallelized
  computations.

- progress_bar:

  should a progress bar be shown. Default is `progress_bar = TRUE`.

- sep:

  a file separator that is used to store maps as csv file.

- path:

  a character string specifying the path to the corresponding orthology
  inference tool. For "BH" and "RBH": path to BLAST.

- ...:

  additional parameters that shall be passed to
  [`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md).

## Details

Given a query organism and a set of subject organsisms that are stored
in the same folder, this function crawls through all subject organsism
and infers the lncRNA homologs in pairwise species comparisons.

## Note

According to Sarropoulos, I., et al. (2019) orthology detection of
lncRNAs was performed by reciprocal BLAST searches. Significant hits
with an e-value \<= 10-3 were selected having an alignment identity \>=
10% OR a minimum alignment length \>= 50 nucleotides.

## References

Sarropoulos I, Marin R, Cardoso-Moreira M, Kaessmann H (2019).
“Developmental dynamics of lncRNAs across mammalian organs and species.”
*Nature*, **571**, 510–514.

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
# example using classic blastn searches
map_generator_lnc(
   query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
   subjects_folder = system.file('seqs/map_gen_example', package = 'orthologr'),
   output_folder   = "orthologs_lncrna",
   comp_cores      = 1
)
# example using  discontiguous megablast used to find somewhat distant sequences
map_generator_lnc(
   query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
   subjects_folder = system.file('seqs/map_gen_example', package = 'orthologr'),
   output_folder   = "orthologs_lncrna",
   task            = "dc-megablast",
   comp_cores      = 1
)
} # }

if (FALSE) { # \dontrun{
# parameter settings based on Sarropoulos, I., et al. (2019)
map_generator_lnc(
   query_file,
   subjects_folder,
   eval                  = 1E-3,
   ortho_detection       = "RBH",
   output_folder,
   min_qry_coverage_hsp  = 0,
   min_qry_perc_identity = 10,
   logical_connective    = "OR",
   min_alig_length       = 50)
} # }
```
