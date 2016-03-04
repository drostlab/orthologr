## ----eval=FALSE----------------------------------------------------------
#  # in case the default execution path of clustalw runs properly on your system
#  multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#            tool    = "clustalw",
#            get_aln = TRUE)
#  

## ----eval=FALSE----------------------------------------------------------
#  # running clustalw using additional parameters
#  # details: system("clustalw2 -help")
#  multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#            tool    = "clustalw",
#            get_aln = TRUE,
#            params  = "-PWMATRIX=BLOSUM -TYPE=PROTEIN")
#  

## ----eval=FALSE----------------------------------------------------------
#  # in case the default execution path of muscle runs properly on your system
#  multi_aln(file    = system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#            tool    = "muscle",
#            get_aln = TRUE)
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  multi_aln(system.file('seqs/aa_seqs.fasta', package = 'orthologr'),
#            tool = "clustalw", get_aln = TRUE)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  multi_aln(system.file('seqs/multi_aln_example.fasta', package = 'orthologr'),
#            tool = "clustalw", get_aln = TRUE)
#  

