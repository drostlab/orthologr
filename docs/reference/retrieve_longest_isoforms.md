# Retrieve the longest isoforms from a proteome file and save results as fasta file

Based on a fasta file storing the peptide isoforms of gene loci and an
annotation file in `gtf` file format, this function extracts the longest
isoform per gene locus and stores the results in a new `fasta` file.
This procedure enables easier downstream analyses such as orthology
inference etc when dealing with proteome `fasta` files which usually
include isoform peptides.

## Usage

``` r
retrieve_longest_isoforms(
  proteome_file,
  annotation_file,
  new_file,
  annotation_format = "gff"
)
```

## Arguments

- proteome_file:

  file path to proteome in `fasta` file format.

- annotation_file:

  file path to the corresponding annotation file in `gtf` file format.

- new_file:

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
# retrieve example data from ENSEMBLGENOMES
proteome <- biomartr::getProteome(db = "refseq", organism = "Arabidopsis thaliana")
annotation <- biomartr::getGFF(db = "refseq", organism = "Arabidopsis thaliana")
# retrieve longest isoforms and store in new file
retrieve_longest_isoforms(proteome_file = proteome, 
                          annotation_file = annotation, 
                          new_file = "Athaliana_pep_longest.fa")
# import new file into R session                          
Athaliana_pep_longest <- Biostrings::readAAStringSet("Athaliana_pep_longest.fa")
} # }
```
