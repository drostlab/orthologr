#!/usr/local/perl

#### Licence #####
# This script is freely available under GNU GPL v3 
# Licence and included in the KaKs_Calculcator 
# project that can be found at 
# https://code.google.com/p/kaks-calculator/
################## 

if(@ARGV<1) {
	die ("Usage: [Fasta file with aligned pairwise sequences]\n");
}
print "*****************************************************************
Function: Parse fasta file with aligned pairwise sequences into AXT file
Reference: Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator: Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263. 
Web Link: Documentation, example and updates at <http://code.google.com/p/kaks-calculator>
*****************************************************************
\n";

print "Input file = $ARGV[0]\n";
print "Opening the input file...\n";
unless (open (MYFILE, $ARGV[0])) {
    die ("Cannot open file\n");
}
@array = <MYFILE>;
chomp @array;

$count = 0;
@seq = ();
@name = ();
$tmp_seq = "";

print "Reading sequences";
while ($count < @array) {
	print ".";
	if(index($array[$count], ">")==0) {
		push(@name, substr($array[$count], 1, length($array[$count])-1));
   		if ($tmp_seq ne "") {
			push(@seq, $tmp_seq);
			$tmp_seq = "";
		}	
	}
	else {
		$tmp_seq .= $array[$count];
	}
	$count++;

	if ($count==@array) {
           push(@seq, $tmp_seq);
        }

}

print "\nFormatting AXT sequences...\n";
$output  = join("-", @name)."\n";
$output .= join("\n", @seq)."\n";
$output .= "\n";

$outfile = $ARGV[1].".axt";
print "Outputing AXT file: $outfile\n";
open(OUTFILE, ">$outfile");
print OUTFILE ($output);

print "\nMission accomplished.\n";
