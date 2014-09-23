/*********************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
*
* Filename: NG86.cpp
* Abstract: Definition of NG86 class.
*
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.21, 2005

  Reference: Nei M, Gojobori T  (1986)  Simple methods for
  estimating the numbers of synonymous and nonsynonymous 
  nucleotide substitutions. Mol Biol Evol 3:418-426.
**********************************************************/

#include "NG86.h"

NG86::NG86() {
	name = "NG";	
}


//Count synonymous(S) and nonsynonymous(N) sites
void NG86::getCondonSite(string codon) {

	int i,j,stop;
	string temp="";
	double syn = 0.0;
	
	if (getAminoAcid(codon)=='!')
		return;
	
	/* Synonymous sites only occur at first and third position in a codon */
	for(i=0, stop=0; i<3; i+=2) {
		for(j=0; j<4; j++) {
			temp = codon;
			if (j!=convertChar(temp[i])) {
				temp[i] = convertInt(j);
				if (getAminoAcid(temp)=='!') {
					stop++;
				}
				else {
					if (getAminoAcid(temp)==getAminoAcid(codon))
						syn++;
				}
			}
		}
	}
	S += (syn/3.0);
	N += (3-stop/3.0-syn/3.0);
}

//Count synonymous(Sd) and nonsynonymous(Nd) differences
void NG86::getCondonDifference(string codon1, string codon2) {

	int i,j,k,diff[CODONLENGTH];
	int num = 0;
	int stop = 0;
	int path = 1;
	double sd_temp = 0.0;
	double nd_temp = 0.0;
	string temp1, temp2;


	if (getAminoAcid(codon1)=='!' || getAminoAcid(codon2)=='!')
		return;

	for (i=0; i<CODONLENGTH; i++) {
		diff[i] = -1;
		if (codon1[i]!=codon2[i])
			diff[num++] = i;
	}

	//two codons are same
	if (num==0) return;
	
	snp += num;

	//Pathway of evolution from the differences
	for (i=1; i<=num; i++) 
		path *= i;

	//Only one difference between two codons
	if (num==1) {
		if (getAminoAcid(codon1)==getAminoAcid(codon2))
			sd_temp++;
		else
			nd_temp++;
	}

	//Two differences between two codons
	if (num==2) {
		for(i=0;i<num;i++)
			for(j=0;j<num;j++)
				if(i!=j) {
					temp1 = codon1;
					temp1[diff[i]] = codon2[diff[i]];
					if (getAminoAcid(temp1)!='!') {
						
						//codon1<->temp1
						if (getAminoAcid(temp1)==getAminoAcid(codon1)) sd_temp++;
						else nd_temp++;

						//temp1<->codon2
						if (getAminoAcid(temp1)==getAminoAcid(codon2)) sd_temp++;
						else nd_temp++;
					}
					else
						stop++;
				}
	}
	//Three differences between two codons
	if(num==3) {
		for (i=0;i<3;i++) {
        	for (j=0;j<3;j++) {
        		for (k=0;k<3;k++) {    				
					if ((i!=j) && (i!=k) && (j!=k)) {						
						temp1 = codon1;
						temp1[diff[i]] = codon2[diff[i]];
						temp2 = temp1;
						temp2[diff[j]] = codon2[diff[j]];
						if (getAminoAcid(temp1)!='!' && getAminoAcid(temp2)!='!') {
							//codon1<->temp1
							if (getAminoAcid(temp1)==getAminoAcid(codon1)) sd_temp++;
							else nd_temp++;

							//temp1<->temp2
							if (getAminoAcid(temp2)==getAminoAcid(temp1)) sd_temp++;
							else nd_temp++;

							//temp2<->codon2
							if (getAminoAcid(codon2)==getAminoAcid(temp2)) sd_temp++;
							else nd_temp++;	

						}					
						else
							stop++;
					}
				}
			}
		}
	}
	if (path==stop) {
		//All pathways are through stop codons
		if (num==2) {
			Sd+=0.5;	Nd+=1.5;
		}
		else {
			Sd+=1.0;	Nd+=2.0; 
		}
	}
	else {
		Sd += (double)(sd_temp/(path-stop));
		Nd += (double)(nd_temp/(path-stop));
	}
	
}

void NG86::PreProcess(string seq1, string seq2) {

	long i;
	
	//Count sites and differences
	for(i=0; i<seq1.length(); i=i+3) {
		getCondonSite(seq1.substr(i,3));
		getCondonSite(seq2.substr(i,3));
		getCondonDifference(seq1.substr(i,3), seq2.substr(i,3));
	}
	
	S/=2.0;
	N/=2.0;

	//Scale the sum of S+N to the length of sequence.
	double y=seq1.length()/(S+N); 
	S*=y;
	N*=y;
}

/* Jukes & Cantor's one-parameter formula for correction */
double NG86::kaks_formula(double p) {
	//Equation (3) in the reference of NG86
	double d = 1-(4*p)/3;
	if (d<0.0) {
		d = NA;
	}
	else {
		d = log(d);
		if (d>0.0)
			d = NA;
		else
			d = (-3.0)*d/4.0;
	}

	return d;
}

string NG86::Run(string seq1, string seq2) {
	
	PreProcess(seq1, seq2);

	Ks = kaks_formula(Sd/S);
	Ka = kaks_formula(Nd/N);

	t = (S*Ks+N*Ka)/(S+N);

	return 	parseOutput();
}



/***********************************************************
  NONE: an in-house algorithm in BGI for testing Ka and Ks.
  NONE is NG86 without correction for multiple substitutions.
************************************************************/
NONE::NONE() {
	name = "NONE";	
}

string NONE::Run(string seq1, string seq2) {
	
	PreProcess(seq1, seq2);

	Ks = Sd/S;
	Ka = Nd/N;

	t = (S*Ks+N*Ka)/(S+N);

	return 	parseOutput();
}
