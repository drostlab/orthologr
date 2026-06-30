# Internal function for dNdS computations

This function takes a pairwise alignment as input file and estimates the
dNdS ratio of the corresponding alignment. Nevertheless, this function
is a helper function for
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md). For
dNdS computations you should use the function:
[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md).

## Usage

``` r
substitutionrate(
  file,
  est.method,
  format = "fasta",
  quiet = FALSE,
  kaks_calc.params = NULL,
  kaks_calc_path = NULL,
  subst_name = NULL
)
```

## Arguments

- file:

  a character string specifying the path to a codon alignment file

- est.method:

  a character string specifying the dNdS estimation method, e.g.
  "Comeron","Li" . Note, that when using "Comeron" as dNdS estimation
  method, the program 'gestimator' is used to compute the corresponding
  dNdS values from a given alignment. The program 'gestimator' can only
  read "fasta" files, hence it is important to use format = "fasta" when
  choosing est.method = "Comeron".

- format:

  a character string specifying the file format in which the alignment
  is stored: "mase", "clustal", "phylip", "fasta" , "msf"

- quiet:

  a logical value specifying whether the output of the coresponding
  interface shall be printed out.

- kaks_calc.params:

  a character string storing additional parameters for KaKs_Claculator
  1.2 . Default is `NULL`. Example: `kaks_calc.params` = "-m NG -m YN".

- kaks_calc_path:

  a character string specifying the execution path to KaKs_Calculator.
  Default is `kaks_calc_path` = `NULL` (meaning that KaKs_Calculator is
  stored and executable in your default `PATH`).

- subst_name:

  a character string specifying the substitution name that shall be
  added to the internal folder path naming. Default is `subst_name` =
  `NULL`.

## Value

A data.table storing the query_id, subject_id, dN, dS, and dNdS values
or a data.table storing the query_id, method, dN, dS, and dNdS values
when using KaKs_Calculator. If the dNdS value cannot be calculated NA is
returned. This can happen because of constraints of the used model. As
each program throws different exception values we set all of them to NA
instead.

## Details

This function takes a pairwise alignments file as input and estimates
dNdS ratios of the corresponding alignments using predefined dNdS
estimation methods.

The dNdS estimation methods available in this function are:

- "Li" : Li's method (1993) -\> provided by the ape package

- "Comeron" : Comeron's method (1995)

dNdS estimation methods provided by KaKs_Calculator 1.2 :

Approximate Methods:

- "NG": Nei, M. and Gojobori, T. (1986)

- "LWL": Li, W.H., et al. (1985)

- "LPB": Li, W.H. (1993) and Pamilo, P. and Bianchi, N.O. (1993)

- "MLWL": (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al.
  (2004)

- "YN": Yang, Z. and Nielsen, R. (2000)

- "MYN": (Modified YN): Zhang, Z., et al. (2006)

- "GMYN": Wang, D.P., et al. Biology Direct. (2009)

- "GY": Goldman, N. and Yang, Z. (1994)

- "MS": (Model Selection)

- "MA": (Model Averaging): based on a set of candidate models,
  Posada, D. (2003)

Maximum-Likelihood Methods:

- GY: Goldman, N. and Yang, Z. (1994)

## References

Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and
nonsynonymous substitution. J. Mol. Evol., 36:96-99.

Charif, D. and Lobry, J.R. (2007) SeqinR 1.0-2: a contributed package to
the R project for statistical computing devoted to biological sequences
retrieval and analysis.

Thornton, K. (2003) libsequence: a C++ class library for evolutionary
genetic analysis. Bioinformatics 19(17): 2325-2327

Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator:
Calculating Ka and Ks through model selection and model averaging.
Genomics Proteomics Bioinformatics 2006 , 4:259-263.

https://code.google.com/archive/p/kaks-calculator

## See also

[`dNdS`](https://drostlab.github.io/orthologr/reference/dNdS.md),
[`multi_aln`](https://drostlab.github.io/orthologr/reference/multi_aln.md),
[`codon_aln`](https://drostlab.github.io/orthologr/reference/codon_aln.md),
[`blast_best`](https://drostlab.github.io/orthologr/reference/blast_best.md),
[`blast_rec`](https://drostlab.github.io/orthologr/reference/blast_rec.md),
[`read.cds`](https://drostlab.github.io/orthologr/reference/read.cds.md)

## Author

Hajk-Georg Drost and Sarah Scharfenberg

## Examples

``` r
if (FALSE) { # \dontrun{

# estimate the dNdS rate using Li's method
substitutionrate(
   file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
   est.method = "Li", 
   format     = "fasta")
 
# estimate the dNdS rate using Comeron's method
substitutionrate(
   file       = system.file("seqs/pal2nal.aln", package = "orthologr"),
   est.method = "Comeron", 
   format     = "fasta")
                 
# estimate the dNdS rate using model averaging provided by
# the KaKs_Calculator 1.2 program
 substitutionrate(
    file       = system.file("seqs/pal2nal.aln", package = "orthologr"), 
    est.method = "MA", 
    format     = "fasta") 
                  
 # estimate the dNdS rate using Nei and Gojobori's method provided by the 
 # KaKs_Calculator 1.2 program
 substitutionrate(
    file       = system.file("seqs/pal2nal.aln", package = "orthologr"), 
    est.method = "NG", 
    format     = "fasta")     
  
  # estimate the dNdS rate using Nei and Gojobori's method AND Yang and Nielsen's 
  # method provided by the KaKs_Calculator 1.2 program 
  # for this purpose we choose: 
  # est.method = "kaks_calc" and kaks_calc.params = "-m NG -m YN"                 
 substitutionrate(
   file             = system.file("seqs/pal2nal.aln", package = "orthologr"),
   est.method       = "kaks_calc", 
   format           = "fasta",
   kaks_calc.params = "-m NG -m YN")
                  
                                          
} # }
```
