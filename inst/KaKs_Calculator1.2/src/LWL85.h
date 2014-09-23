/************************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: LWL85.h
* Abstract: Declaration of LWL85 and Modified LWL85 class.

* Version: 1.0
* Author: Zhang Zhang(zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005

  References: 
  Li WH, Wu CI, Luo CC  (1985)  A new method for
  estimating synonymous and nonsynonymous rates of nucleotide 
  substitution considering the relative likelihood of nucleotide
  and codon changes. Mol. Biol. Evol. 2:150-174.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
*************************************************************/

#if !defined(LWL85_H)
#define  LWL85_H

#include "base.h"

using namespace std;

/* Nondegenerate, twofold and fourfold when considering position 1,2,3 of 64 Codon */
/* 0: 0-fold,	2: 2-fold,		4: 4-fold */

/* The following is only for standard genetic code.
   Considering more genetic codes, we use the function of "getCodonClass".
const int CodonClassDict[]={0,0,2,  0,0,2,	2,0,2,  2,0,2,	0,0,4,  0,0,4,	0,0,4,  0,0,4,
							0,0,2,  0,0,2,	0,0,0,  0,0,0,	0,0,2,  0,0,2,	0,0,0,  0,0,0,
							0,0,4,  0,0,4,	2,0,4,  2,0,4,	0,0,4,  0,0,4,	0,0,4,  0,0,4,
							0,0,2,  0,0,2,	0,0,2,  0,0,2,	0,0,4,  0,0,4,	2,0,4,  2,0,4,
							0,0,2,  0,0,2,	0,0,2,  0,0,0,	0,0,4,  0,0,4,	0,0,4,  0,0,4,
							0,0,2,  0,0,2,	0,0,2,  0,0,2,	0,0,2,  0,0,2,	2,0,2,  2,0,2,
							0,0,4,  0,0,4,	0,0,4,  0,0,4,	0,0,4,  0,0,4,	0,0,4,  0,0,4,
							0,0,2,  0,0,2,	0,0,2,  0,0,2,	0,0,4,  0,0,4,	0,0,4,  0,0,4}; */



/* LWL85 class */
class LWL85: public Base {
	
public:
	LWL85();

	/* Main function for calculating Ka&Ks */
	string Run(string seq1, string seq2);

protected:
	/* preprocess in main function of Run */
	void preProcess(string seq1, string seq2);
	/* Calculate synonymous and nonsynonymous sites and differences on two compared codons */
	void CountSiteAndDiff(string str1, string str2);
	/* Return 0,2,or 4 of the codon at a given position */
	int  getCodonClass(string codon, int pos);
	/* Calculate synonymous and nonsynonymous differences of two codons at a given position */
	virtual int  TransitionTransversion(string codon1, string codon2, int pos);

public:

	/* Proportions of transitional(Pi) and transversional(Qi) difference: Pi=Si/Li,  Qi=Vi/Li  */
	double Pi[5], Qi[5];
	/* Numbers of transitional(A) and transversional(B) substitutions per i-th type site */
	double A[5], B[5];

	
protected:
	/* Number of synonymous and nonsynonymous differences of pair-codon */
	double Si_temp[5], Vi_temp[5];	
}; 


class MLWL85: public LWL85 {
	
public:
	MLWL85();

	/* Main function for calculating Ka&Ks */
	string Run(string str1, string str2);
	
protected:
	/* Calculate the transition & transversion between two codons at a given position*/
	int TransitionTransversion(string codon1, string codon2, int pos);
};


#endif


