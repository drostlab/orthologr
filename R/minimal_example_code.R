
# example to check tcoffee run -> much output that we should try to eliminate
aln <- multi_aln(file = "data/aa_seqs.fasta", tool = "tcoffee", get_aln=TRUE)

#example to check pal2nal run
codon_aln <- codon_aln(file_aln="data/aa_seqs.aln", file_nuc = "data/dna_seqs.fasta", 
                       format = "clustal", tool="pal2nal", get_aln=TRUE,
                       path="exec/pal2nal.v14/")