## ---- echo = FALSE, message = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  library(orthologr)
#  
#  # get a dNdS table using:
#  # 1) reciprocal best hit for orthology inference (RBH)
#  # 2) clustalw for pairwise amino acid alignments
#  # 3) pal2nal for codon alignments
#  # 4) Yang, Z. and Nielsen, R. (2000) (YN) for dNdS estimation
#  # 5) single core processing 'comp_cores = 1'
#  dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#        subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#        ortho_detection = "RBH",
#        aa_aln_type     = "multiple",
#        aa_aln_tool     = "clustalw",
#        codon_aln_tool  = "pal2nal",
#        dnds_est.method = "YN",
#        comp_cores      = 1,
#        clean_folders   = TRUE,
#        quiet           = TRUE )
#  
#  

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  library(orthologr)
#  
#  # get dNdS estimated for orthologous genes between A. thaliana and A. lyrata
#  
#  Ath_Aly_dnds <- dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#                        subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#                        dnds_est.method = "YN",
#                        comp_cores      = 1,
#                        clean_folders   = TRUE,
#                        quiet           = TRUE )
#  
#  
#  # filter for:
#  # 1) all dN values having an NA value are omitted
#  # 2) all dS values having an NA value are omitted
#  # 3) all dNdS values >= 2 are omitted
#  
#  filter_dNdS(Ath_Aly_dnds, dnds.threshold = 2)
#  
#  

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  library(orthologr)
#  
#  # get a dNdS table using:
#  # 1) reciprocal best hit for orthology inference (RBH)
#  # 2) pairwise amino acid alignments using Needleman-Wunsch
#  # 3) pal2nal for codon alignments
#  # 4) Comeron (1995) for dNdS estimation
#  # 5) single core processing 'comp_cores = 1'
#  dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#        subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#        ortho_detection = "RBH",
#        aa_aln_type     = "pairwise",
#        aa_aln_tool     = "NW",
#        codon_aln_tool  = "pal2nal",
#        dnds_est.method = "Comeron",
#        comp_cores      = 1,
#        clean_folders   = TRUE,
#        quiet           = TRUE )
#  
#  

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  library(orthologr)
#  
#  # using the `aa_aln_path` or `blast_path` arguments
#  dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#        subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#        ortho_detection = "RBH",
#        blast_path      = "here/path/to/blastp",
#        aa_aln_type     = "multiple",
#        aa_aln_tool     = "clustalw",
#        aa_aln_path     = "here/path/to/clustalw",
#        codon_aln_tool  = "pal2nal",
#        dnds_est.method = "Comeron",
#        comp_cores      = 1,
#        clean_folders   = TRUE,
#        quiet           = TRUE )
#  
#  

## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  library(orthologr)
#  
#  # get dNdS estimated for orthologous genes between A. thaliana and A. lyrata
#  # using additional parameters:
#  
#  # get a dNdS table using:
#  # 1) reciprocal best hit for orthology inference (RBH)
#  # 2) multiple amino acid alignments using MAFFT
#  # 3) pal2nal for codon alignments
#  # 4) Comeron (1995) for dNdS estimation
#  # 5) single core processing 'comp_cores = 1'
#  Ath_Aly_dnds <- dNdS( query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#                        subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#                        ortho_detection = "RBH",
#                        blast_params    = "-matrix BLOSUM80",
#                        aa_aln_tool     = "mafft",
#                        aa_aln_params   = "--maxiterate 1 --clustalout",
#                        dnds_est.method = "Comeron",
#                        comp_cores      = 1,
#                        clean_folders   = TRUE,
#                        quiet           = TRUE )
#  
#  
#  # filter for:
#  # 1) all dN values having an NA value are omitted
#  # 2) all dS values having an NA value are omitted
#  # 3) all dNdS values >= 0.1 are omitted
#  
#  filter_dNdS(Ath_Aly_dnds, dnds.threshold = 0.1)
#  
#  

