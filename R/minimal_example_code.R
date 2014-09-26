
# example to check tcoffee run -> much output that we should try to eliminate
#aln <- multi_aln(file = "data/aa_seqs.fasta", tool = "tcoffee", get_aln = TRUE, path = "")

#example to check pal2nal run -> use your own path to pal2nal
# codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#                        file_nuc = system.file('seqs/data/dna_seqs.fasta', package = 'orthologr'), 
#                        format = "clustal", tool = "pal2nal", get_aln = TRUE)


# dNdS(query_file = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#      subject_file = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#      ortho_detection = "RBH", blast_path = NULL,
#      multialn_tool = "t_coffee", multialn_path = "",
#      codonaln_tool = "pal2nal", dnds_est.method = "YN", comp_cores = 2, quiet = FALSE)
# # 
# 
# compute_dnds <- function(x, 
#                          multialn_tool="clustalw", multialn_path = NULL,
#                          codonaln_tool="pal2nal", dnds_est.method = "Comeron",
#                          quiet = FALSE)