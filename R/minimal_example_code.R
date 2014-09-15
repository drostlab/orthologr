
# example to check tcoffee run -> much output that we should try to eliminate
aln <- multi_aln(file = "data/aa_seqs.fasta", tool = "tcoffee", get.aln=TRUE)

#example to check pal2nal run, change path!
codon_aln <- codon_aln(file.aln="data/aa_seqs.aln", file.nuc = "data/dna_seqs.fasta", 
                       format = "clustal", tool="pal2nal", get.aln=TRUE,
                       path="/home/sarah/Programs/pal2nal.v14/")