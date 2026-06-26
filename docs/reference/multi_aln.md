# Compute Multiple Sequence Alignments

This function takes a FASTA file containing DNA or amino acid sequences
that shall be aligned and computes a multiple alignment using a defined
multiple alignment tool.

## Usage

``` r
multi_aln(
  file,
  tool,
  get_aln = FALSE,
  path = NULL,
  multi_aln_name = NULL,
  params = NULL,
  quiet = FALSE,
  clean_folders = FALSE
)
```

## Arguments

- file:

  a character string specifying the path to the file storing the
  sequences in FASTA format.

- tool:

  a character string specifying the program that should be used:
  "clustalw", "t_coffee", "muscle", "clustalo", and "mafft".

- get_aln:

  a logical value indicating whether the produced alignment should be
  returned.

- path:

  a character string specifying the path to the multiple alignment
  program (in case you don't use the default path).

- multi_aln_name:

  a character string specifying the name of the stored alignment file.
  Default is `multi_aln_name` = `NULL` denoting a default name:
  'toolname.aln' .

- params:

  a character string listing the input paramters that shall be passed to
  corresponding alignment tool specified in the `tool` argument. Default
  is `NULL`, implicating that a set of default parameters is used when
  running the correponding tool, e.g. clustalw. Example: `params` =
  "-PWMATRIX=BLOSUM -TYPE=PROTEIN" for `tool` = "clustalw.

- quiet:

  a logical value specifying whether the output of the corresponding
  multiple alignment tool shall be printed out to the console. Default
  is `quiet` = `FALSE`.

- clean_folders:

  a boolean value spefiying whether all internall folders storing the
  output of used programs shall be removed. Default is `clean_folders` =
  `FALSE`.

## Value

In case the argument `get_aln` is set `TRUE`, an object of class
alignment of the seqinr package is returned.

## Details

This function provides an interface between R and common multiple
alignment programs such as "clustalw", "t_coffee", "muscle", "clustalo",
and "mafft".

- CLUSTALW :

  Different operating systems perform different execution calls to the
  clustalw program:

  - MacOS: 'clustalw2'

  - Linux: 'clustalw'

  - Windows: 'clustalw2.exe'

  In case you use the default path to the clustalw program, depending on
  your operating system, the following calls to clustalw should work
  properly on your system:

  - MacOS: `system("clustalw2 -help")`

  - Linux: `system("clustalw -help")`

  - Windows: `system("clustalw2.exe -help")`

  In case these procedures don't work properly, please use the `path`
  argument to specify the 'clustalw' execution path on your system:

  - MacOS: `system("path/to/clustalw/clustalw2 -help")`

  - Linux: `system("path/to/clustalw/clustalw -help")`

  - Windows: `system("path/to/clustalw/clustalw2.exe -help")`

- T_COFFEE :

  In case you use the default path to the t_coffee program, the
  following calls to clustalw should work properly on your system:

  `system("t_coffee -version")`

  In case this procedures doesn't work properly, please use the `path`
  argument to specify the 't_coffee' execution path on your system:

  `system("path/to/t_coffee/t_coffee -version")`

- MUSCLE :

  In case you use the default path to the muscle program, the following
  calls to muscle should work properly on your system:

  `system("muscle -help")`

  In case this procedures doesn't work properly, please use the `path`
  argument to specify the 'muscle' execution path on your system:

  `system("path/to/muscle/muscle -help")`

- CLUSTALO :

  In case you use the default path to the clustalo program, the
  following calls to clustalo should work properly on your system:

  `system("clustalo --help")`

  In case this procedures doesn't work properly, please use the `path`
  argument to specify the 'clustalo' execution path on your system:

  `system("path/to/clustalo/clustalo --help")`

- MAFFT :

  In case you use the default path to the mafft program, the following
  calls to mafft should work properly on your system:

  `system("mafft -help")`

  In case this procedures doesn't work properly, please use the `path`
  argument to specify the 'mafft' execution path on your system:

  `system("path/to/mafft/mafft -help")`

## Note

Note that when using the `clustalw.params`, ... parameters, make sure
the corresponding alignment tool returns a file in clustal format
\*.aln. This is only important when `get_aln` = `TRUE`.

## References

- CLUSTALW:

  Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA,
  McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD,
  Gibson TJ, Higgins DG. (2007). Clustal W and Clustal X version 2.0.
  Bioinformatics, 23, 2947-2948.

  <https://www.clustal.org/clustal2/>

  <https://www.ebi.ac.uk/jdispatcher/>

- T_COFFEE

  T-Coffee: A novel method for multiple sequence alignments. Notredame,
  Higgins, Heringa, JMB, 302(205-217). 2000.

  <https://tcoffee.readthedocs.io/en/latest/>

  <https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html>

- MUSCLE:

  Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high
  accuracy and high throughput. Nucleic Acids Res. 32(5):1792-1797.

  Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with
  reduced time and space complexity. BMC Bioinformatics, (5) 113.

  <https://www.drive5.com/muscle/>

  <https://www.drive5.com/muscle/manual/>

- CLUSTALO:

  Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R,
  McWilliam H, Remmert M, Soeding J, Thompson JD, Higgins DG (2011).
  Fast, scalable generation of high-quality protein multiple sequence
  alignments using Clustal Omega. Molecular Systems Biology 7:539
  doi:10.1038/msb.2011.75

  <https://www.ebi.ac.uk/jdispatcher/msa/clustalo/>

  <https://www.genome.jp/tools-bin/clustalw>

- MAFFT :

  Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)
  MAFFT multiple sequence alignment software version 7: improvements in
  performance and usability.

  <https://mafft.cbrc.jp/alignment/software/>

  <https://mafft.cbrc.jp/alignment/software/manual/manual.html>

  <https://mafft.cbrc.jp/alignment/software/tips0.html>

## Author

Hajk-Georg Drost and Sarah Scharfenberg

## Examples

``` r
if (FALSE) { # \dontrun{

# CLUSTALW Example:

# in case the default execution path of clustalw runs properly on your system
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "clustalw", 
          get_aln = TRUE)

# in case the default execution path of clustalw is not set within the default path
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'), 
          tool    = "clustalw", 
          get_aln = TRUE, 
          path    = "path/to/clustalw/")

# running clustalw using additional parameters
# details: system("clustalw2 -help")
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "clustalw", 
          get_aln = TRUE, 
          params  = "-PWMATRIX=BLOSUM -TYPE=PROTEIN")
          
          
# T_COFFEE Example:

# in case the default execution path of t_coffee runs properly on your system
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "t_coffee", 
          get_aln = TRUE)

# in case the default execution path of t_coffee is not set within the default path
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "t_coffee", 
          get_aln = TRUE, 
          path    = "path/to/t_coffee/")

# running t_coffee using additional parameters
# details: https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "t_coffee", 
          get_aln = TRUE,
          params  = "-mode expresso")
          


# MUSCLE Example:  

# in case the default execution path of muscle runs properly on your system
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "muscle", 
          get_aln = TRUE)

# in case the default execution path of muscle is not set within the default path
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "muscle", 
          get_aln = TRUE, 
          path    = "path/to/muscle/")

# running muscle using additional parameters
# details: https://www.drive5.com/muscle/manual/
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "muscle", 
          get_aln = TRUE,
          params  = "-diags -clwstrict") 
          

# CLUSTALO Example:  

# in case the default execution path of clustalo runs properly on your system
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "clustalo", 
          get_aln = TRUE)

# in case the default execution path of clustalo is not set within the default path
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "clustalo", 
          get_aln = TRUE, 
          path    = "path/to/clustalo/")

# running clustalo using additional parameters
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "clustalo", 
          get_aln = TRUE,
          params  = "--outfmt clu")         
          
          
                                            
# MAFFT Example:  

# in case the default execution path of mafft runs properly on your system
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "mafft", 
          get_aln = TRUE)

# in case the default execution path of mafft is not set within the default path
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "mafft", 
          get_aln = TRUE, 
          path    = "path/to/mafft/")

# running mafft using additional parameters
# details: https://mafft.cbrc.jp/alignment/software/manual/manual.html
multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
          tool    = "mafft", 
          get_aln = TRUE,
          params  = "--maxiterate 1 --clustalout")         
          
                                                                                  
} # }
```
