# Compute a Codon Alignment

This function takes a protein alignment and the corresponding coding
sequences as fasta files and computes the corresponding codon alignment
based on the models specified in the PAL2NAL program.

## Usage

``` r
codon_aln(
  file_aln,
  file_nuc,
  format = "clustal",
  tool = "pal2nal",
  params = NULL,
  codon_aln_name = NULL,
  get_aln = FALSE,
  quiet = FALSE
)
```

## Arguments

- file_aln:

  a character string specifying the path to the file storing the protein
  alignment in CLUSTAL or FASTA format.

- file_nuc:

  a character string specifying the path to the file storing the coding
  sequences in multiple FASTA format.

- format:

  a character string specifying the file format used to store the codon
  alignment, e.g. "fasta", "clustal".

- tool:

  a character string specifying the program that should be used e.g.
  "pal2nal".

- params:

  a character string specifying additional parameters that shall be
  passed to PAL2NAL. For example `params` = `"-codontable 2"` would
  specify a codon table: *vertebrate mitochondrial code* instead if the
  *universal code* (default). The default value is `params` = `NULL`.
  Default is `codon_aln_name` = `NULL` denoting a default name:
  'tool_name.aln' .

- codon_aln_name:

  a character string specifying the name of the stored alignment file.

- get_aln:

  a logical value indicating whether the produced alignment should be
  returned.

- quiet:

  a logical value specifying whether a successful interface call shall
  be printed out.

## Value

In case `get_aln` = `TRUE` an seqinr alignment object is returned.

## Details

This function provides an interface between the R language and common
codon alignment tools such as "PAL2NAL". Codon alignments can be used to
quantify the evolutionary pressure acting on protein sequences based on
dNdS estimation.

The PAL2NAL program is used to convert a multiple sequence alignment or
pairwise alignment of proteins and the corresponding DNA (or mRNA)
sequences into a codon-based DNA alignment. The output codon-based DNA
alignment can then be used for
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md)
estimation.

The popular codon alignment tool PAL2NAL is included in this package and
is being called when choosing `tool` = `"pal2nal"`.

## Note

The PAL2NAL program automatically checks

- If you use the same IDs are used in the protein alignment file and DNA
  (or mRNA) file -\> in this case the sequences don't have to be in the
  same order

- If you don't use the same IDs are used in the protein alignment file
  and DNA (or mRNA) file -\> you have to rearrange the IDs or sequences
  so that sequences in the protein alignment file and DNA (or mRNA) file
  have the same order

## References

Mikita Suyama, David Torrents, and Peer Bork (2006) PAL2NAL: robust
conversion of protein sequence alignments into the corresponding codon
alignments. Nucleic Acids Res. 34, W609-W612.

https://www.bork.embl.de/services.html

https://github.com/abacus-gene/paml/wiki

## See also

[`pairwise_aln`](https://drostlab.github.io/orthologr/reference/pairwise_aln.md),
[`multi_aln`](https://drostlab.github.io/orthologr/reference/multi_aln.md),
[`substitutionrate`](https://drostlab.github.io/orthologr/reference/substitutionrate.md),
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md),
[`divergence_stratigraphy`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)

## Author

Sarah Scharfenberg and Hajk-Georg Drost

## Examples

``` r

if (FALSE) { # \dontrun{

# performing a codon alignment using PAL2NAL
codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
                       file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'), 
                       format   = "clustal", 
                       tool     = "pal2nal", 
                       get_aln  = TRUE)
                       
                       
# running PAL2NAL with additional parameters (e.g. a different codon table)
codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
                       file_nuc = system.file('seqs/dna_seqs.fasta', package = 'orthologr'), 
                       format   = "clustal", 
                       tool     = "pal2nal", 
                       get_aln  = TRUE,
                       params   = "-codontable 2")
                        
} # }
                       
```
