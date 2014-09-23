/*********************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
*
* Filename: LPB93.cpp
* Abstract: Definition of LPB93 and MLPB93 class.
*
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005

  References: 
  Li WH  (1993)  Unbiased estimation of the Rates of synonymous
  and nonsynonymous substitution. J. Mol. Evol. 36:96-99.

  Pamilo P, Bianchi NO  (1993)  Evolution of the Zfx and Zfy 
  genes: rates and interdependence between the genes. Mol. Biol.
  Evol. 10:271-281.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
**********************************************************/


#include "LPB93.h"

LPB93::LPB93() {
	name = "LPB";
}

/* Similar to LWL85 except the formulas for calculating ka and ks*/
string LPB93::Run(string seq1, string seq2) {

	preProcess(seq1, seq2);
	
	Ks = B[4] + (L[2]*A[2] + L[4]*A[4])/(L[2] + L[4]);
	Ka = A[0] + (L[0]*B[0] + L[2]*B[2])/(L[0] + L[2]);
	
	Sd = L[2]*A[2] + L[4]*K[4];
	Nd = L[2]*B[2] + L[0]*K[0];
	
	S = Sd/Ks;
	N = Nd/Ka;

	t = (L[0]*K[0]+L[2]*K[2]+L[4]*K[4])/(L[0]+L[2]+L[4]);

	return parseOutput();
}


/*The difference between LPB93 and MLPB93 focuses on the definition of transition & transversion*/
MLPB93::MLPB93() {
	name = "MLPB";
}

/*For more detail see reference: Tzeng Y-H, Pan R, Li W-H  (2004)  Mol. Biol. Evol.*/
int MLPB93::TransitionTransversion(string codon1, string codon2, int pos) {
	
	//1:synonymous, 0:nonsynonymous, -1:uncalculate
	int isSyn=-1;
	
	//Ile: ATT, ATC, ATA; Met: ATA
	if ((codon1=="ATA" && codon2=="ATG" && pos==2)||(codon1=="ATG" && codon2=="ATA" && pos==2)) {
		isSyn = 0;		
	}
	if ((codon1=="ATA" && (codon2=="ATC"||codon2=="ATT") && pos==2) || ( (codon1=="ATC"||codon1=="ATT") && codon2=="ATA" && pos==2)) {
		isSyn = 1;		
	}
	
	//Arg: CGA, CGG, AGA, AGG
	if ((codon1=="CGA" && codon2=="AGA" && pos==0)||(codon1=="AGA" && codon2=="CGA" && pos==0)) {
		isSyn = 1;
	}
	if ((codon1=="CGG" && codon2=="AGG" && pos==0)||(codon1=="AGG" && codon2=="CGG" && pos==0)) {
		isSyn = 1;
	}

	//Synonymous: A<->G, C<->T
	//Normal situation	
	if (isSyn==-1) {
		int sum = convertChar(codon1[pos]) + convertChar(codon2[pos]);	
		if (sum==5 || sum==1)
			isSyn = 1;
		else
			isSyn = 0;
	}
	
	int class1=getCodonClass(codon1,pos);
	int class2=getCodonClass(codon2,pos);
	if (isSyn==1) {
		Si_temp[class1] += 0.5;
		Si_temp[class2] += 0.5;
	}
	if (isSyn==0) {
		Vi_temp[class1] += 0.5;
		Vi_temp[class2] += 0.5;
	}

	return 0;
}




