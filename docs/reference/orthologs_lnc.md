# Orthology Inference of lncRNAs

This function takes nucleotide or protein sequences for a set of
organisms and performs orthology inference to detect orthologous genes
within the given organisms based on selected orthology inference
programs.

This function takes nucleotide or protein sequences for a set of
organisms and performs orthology inference to detect orthologous genes
within the given organisms based on selected orthology inference
programs.

## Usage

``` r
orthologs_lnc(
  query_file,
  subject_file,
  task = "blastn",
  eval = "1E-5",
  ortho_detection = "RBH",
  max.target.seqs = 10000,
  output.path = getwd(),
  comp_cores = 1,
  path = NULL
)

orthologs_lnc(
  query_file,
  subject_file,
  task = "blastn",
  eval = "1E-5",
  ortho_detection = "RBH",
  max.target.seqs = 10000,
  output.path = getwd(),
  comp_cores = 1,
  path = NULL
)
```

## Arguments

- query_file:

  a character string specifying the path to the sequence file of
  interest (query organism).

- subject_file:

  a character string specifying the paths to the sequence files of
  interest (subject organisms).

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

  a numeric value specifying the E-Value cutoff for BLAST hit detection.

- ortho_detection:

  a character string specifying the orthology inference method that
  shall be performed to detect orthologous genes. Options are:

  - `ortho_detection = "RBH"` (BLAST reciprocal best hit) (Default)

  - `ortho_detection = "BH"` (BLAST best hit)

- max.target.seqs:

  a numeric value specifying the number of aligned sequences to keep.
  Please be aware that max.target.seqs selects best hits based on the
  database entry and not by the best e-value. See details here:
  https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166
  .

- output.path:

  path to which output shall be stored.

- comp_cores:

  a numeric value specifying the number of cores to be used for
  multicore computations.

- path:

  a character string specifying the path to the corresponding orthology
  inference tool. For "BH" and "RBH": path to BLAST, "PO": path to
  ProteinOrtho 5.07, "OrthoMCL": path to OrthoMCL.

## Value

A data.table storing the query_ids of orthologous genes in the first
column, the subject_ids of orthologous genes in the second column and
the amino acid sequences in the third column.

A data.table storing the query_ids of orthologous genes in the first
column, the subject_ids of orthologous genes in the second column and
the amino acid sequences in the third column.

## Details

This function takes sequence files of a query organism and a subject
organism and performs orthology inference using a defined orthology
inference method to dectect orthologous genes.

The following interfaces are implemented in the `orthologs` function:

BLAST based methods:

- BLAST best hit (BH)

- BLAST reciprocal best hit (RBH)

This function takes sequence files of a query organism and a subject
organism and performs orthology inference using a defined orthology
inference method to dectect orthologous genes.

The following interfaces are implemented in the `orthologs` function:

BLAST based methods:

- BLAST best hit (BH)

- BLAST reciprocal best hit (RBH)

## References

BLAST: https://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml

ProteinOrtho: https://www.bioinf.uni-leipzig.de/Software/proteinortho/

BLAST: https://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml

ProteinOrtho: https://www.bioinf.uni-leipzig.de/Software/proteinortho/

## See also

[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)

## Author

Hajk-Georg Drost

## Examples

``` r
if (FALSE) { # \dontrun{

### BLAST Reciprocal Best Hit
# perform orthology inference using BLAST reciprocal best hit
# and fasta sequence files storing protein sequences
orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "RBH")
          
          
### BLAST Best Hit
# perform orthology inference using BLAST best hit
# and fasta sequence files storing protein sequences
orthologs.lnc(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "BH")


# multicore version          
orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "RBH", 
          comp_cores      = 2)          
          
          
          
} # }
if (FALSE) { # \dontrun{

### BLAST Reciprocal Best Hit
# perform orthology inference using BLAST reciprocal best hit
# and fasta sequence files storing protein sequences
orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "RBH")
          
          
### BLAST Best Hit
# perform orthology inference using BLAST best hit
# and fasta sequence files storing protein sequences
orthologs.lnc(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "BH")


# multicore version          
orthologs.lnc(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
          subject_file   = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
          ortho_detection = "RBH", 
          comp_cores      = 2)          
          
          
          
} # }
```
