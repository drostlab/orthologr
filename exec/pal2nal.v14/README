#=================================================================#
#                                                                 #
#  PAL2NAL: robust conversion of protein sequence alignments      #
#           into the corresponding codon-based DNA alignments     #
#                                                                 #
#               version 1.0        June 15, 2005                  #
#               version 1.1        October 11, 2005               #
#               version 2.1        March 28, 2006                 #
#               version 8.0        May 29, 2006                   #
#               version 10.0       August 1, 2006                 #
#               version 11.0       August 9, 2006                 #
#               version 12.0       February 2, 2007               #
#               version 12.1       June 22, 2009                  #
#               version 12.2       September 9, 2009              #
#               version 13.0       July 26, 2010                  #
#               version 14.0       December 2, 2011               #
#                                                                 #
#               Mikita Suyama (mikita@genome.med.kyoto-u.ac.jp)   #
#                                                                 #
#=================================================================#


#------------------#
# What is pal2nal?
#------------------#

PAL2NAL is a program that converts a multiple sequence alignment
of proteins and the corresponding DNA (or mRNA) sequences into
a codon-based DNA alignment. The program automatically assigns
the corresponding codon sequence even if the input DNA sequence
has mismatches with the input protein sequence, or contains UTRs,
polyA tails. It can also deal with frame shifts in the input
alignment, which is suitable for the analysis of pseudogenes.
The resulting codon-based DNA alignment can further be subjected
to the calculation of synonymous (Ks) and non-synonymous (Ka)
substitution rates.


#-----------#
# Reference
#-----------#

If you use PAL2NAL, please cite the following paper:

  - Mikita Suyama, David Torrents, and Peer Bork (2006)
    PAL2NAL: robust conversion of protein sequence alignment into
    the corresponding codon alignments.
    Nucleic Acids Res. 34:W609-W612.


#-------#
# Files
#-------#

The distribution version should contain the following files:

    README                 - This document
    pal2nal.pl             - The script (version 14) (written in Perl)
    test.aln               - test data (protein alignment)
    test.nuc               - test data (DNA sequences)

    for_paml (directory)
      test.cnt             - control file for codeml
      test.tree            - tree file used for codeml
      test.codeml.ori      - an example of codeml output


#-------#
# Usage
#-------#

Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options]

    pep.aln:     protein alignment either in CLUSTAL or FASTA format

                 - works not only a pairwise alignment but also
                   the alignment with more than 2 sequences.

                 - if there are frame shifts in your alignment,
                   you have to specify those positions by numbers:
                   for example '2' in the alignment means that
                   there are only two bases (i.e. one base deletion)
                   (see 'test.aln').


    nuc.fasta:   DNA sequences
                 (single multiple fasta format file,
                  or may be separated files)

    Options:

       -h            Show help

       -output (clustal|paml|fasta|codon)
                     Output format, default = clustal

       -blockonly    Show only user specified blocks
                     '#' under CLUSTAL alignment (see example)

       -nogap        remove columns with gaps and inframe stop codons

       -nomismatch   remove mismatched codons (mismatch between
                     pep and cDNA) from the output

       -codontable (1(default)|2|3|4|5|6|9|10|11|12|13|14|15|16|21|22|23)
                     NCBI GenBank codon table
                     1  Universal code
                     2  Vertebrate mitochondrial code
                     3  Yeast mitochondrial code
                     4  Mold, Protozoan, and Coelenterate Mitochondrial code
                        and Mycoplasma/Spiroplasma code
                     5  Invertebrate mitochondrial
                     6  Ciliate, Dasycladacean and Hexamita nuclear code
                     9  Echinoderm and Flatworm mitochondrial code
                    10  Euplotid nuclear code
                    11  Bacterial, archaeal and plant plastid code
                    12  Alternative yeast nuclear code
                    13  Ascidian mitochondrial code
                    14  Alternative flatworm mitochondrial code
                    15  Blepharisma nuclear code
                    16  Chlorophycean mitochondrial code
                    21  Trematode mitochondrial code
                    22  Scenedesmus obliquus mitochondrial code
                    23  Thraustochytrium mitochondrial code

       -html         HTML output (only for the web server)

       -nostderr     No STDERR messages (only for the web server)



    - The correspondence of IDs between pep.aln and nuc.fasta is
      automatically checked:
        - If you use the same IDs in both pep.aln and nuc.fasta,
          the sequences don't have to be in the same order.
        - If not, the order of the sequences in pep.aln and nuc.fasta
          has to be the same.

    - IDs in pep.aln are used in the output.


Example:  pal2nal.pl  test.aln  test.nuc  -output paml  -nogap


#---------------------------------#
# How to calculate Ks, Ka values?
#---------------------------------#

To calclate Ks, Ka values, you need the codeml program, which
is included in the PAML package. You can download PAML from

  http://abacus.gene.ucl.ac.uk/software/paml.html

As an example, a control file (test.cnt) and a tree file (test.tree)
are in the "for_paml" sub-directory.

These control file and tree file are designed for the 'test' data
used in PAL2NAL.

Example:

   pal2nal.pl  test.aln  test.nuc  -output paml  -nogap  >  for_paml/test.codon

   cd for_paml

   codeml  test.cnt

     You can find the output of codeml in "test.codeml".
     Ks, Ka values are very end of the output file.
     Just for comparison, there is a sample output file, "test.codeml.ori".


#------------#
# WWW server
#------------#

http://www.bork.embl.de/pal2nal
 or
http://www.genome.med.kyoto-u.ac.jp/cgi-bin/suyama/pal2nal/index.cgi


#---------#
# Contact
#---------#

If you have any questions or comments, please email me:

  Mikita Suyama
  Center for Genomic Medicine,
  Graduate School of Medicine, Kyoto University,
  606-8501 Kyoto, JAPAN
  tel: +81 75 753 4383
  fax: +81 75 753 4382
  email: mikita@genome.med.kyoto-u.ac.jp

