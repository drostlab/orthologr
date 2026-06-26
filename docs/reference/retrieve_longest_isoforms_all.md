# Retrieve the longest isoforms of several proteome files stored in a folder

Based on a fasta file storing the peptide isoforms of gene loci and an
annotation file in `gtf` file format, this function extracts the longest
isoform per gene locus and stores the results in a new `fasta` file.
This procedure enables easier downstream analyses such as orthology
inference etc when dealing with proteome `fasta` files which usually
include isoform peptides.

## Usage

``` r
retrieve_longest_isoforms_all(
  proteome_folder,
  annotation_folder,
  output_folder,
  annotation_format = "gff"
)
```

## Arguments

- proteome_folder:

  file path to proteome in `fasta` file format.

- annotation_folder:

  file path to the corresponding annotation file in `gtf` file format.

- output_folder:

  file path to new file storing only peptide sequences of the longest
  isoforms.

- annotation_format:

  format of `annotation_file`. Options are:

  - `annotation_file = "gff"` (default)

  - `annotation_file = "gtf"`

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{
orgs <- c("Arabidopsis lyrata", 
          "Capsella rubella", "Solanum lycopersicum")
# download proteome files for all species          
biomartr::getProteomeSet(db = "refseq", organisms = orgs, path = "of_proteomes")
# download annotation files for all species          
biomartr::getGFFSet(db = "refseq", organisms = orgs, path = "of_gff")
# select longest splice variant per gene locus
retrieve_longest_isoforms_all(proteome_folder = "of_proteomes", 
                              annotation_folder = "of_gff",
                              annotation_format = "gff", 
                              output_folder = "of_proteomes_longest_sv")
} # }
```
