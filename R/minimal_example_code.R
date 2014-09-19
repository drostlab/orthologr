
# example to check tcoffee run -> much output that we should try to eliminate
#aln <- multi_aln(file = "data/aa_seqs.fasta", tool = "tcoffee", get_aln = TRUE, path = "/Users/Hajk/Desktop/Projekte/bioinformatics tools/t_coffee/bin")

#example to check pal2nal run -> use your own path to pal2nal
# codon_aln <- codon_aln(file_aln = system.file('seqs/aa_seqs.aln', package = 'orthologr'),
#                        file_nuc = system.file('seqs/data/dna_seqs.fasta', package = 'orthologr'), 
#                        format = "clustal", tool = "pal2nal", get_aln = TRUE)