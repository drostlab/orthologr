/************************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: LPB93.h
* Abstract: Declaration of LPB93 and Modified LPB93 class.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
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
*************************************************************/

#if !defined(LPB93_H)
#define  LPB93_H

#include "base.h"
#include "LWL85.h"	//derived from LWL85

using namespace std;

/* LPB93 class */
class LPB93: public LWL85 {
	
public:
	LPB93();
	/* Main function of calculating kaks */
	string Run(string seq1, string seq2);
}; 


class MLPB93: public LPB93 {
	
public:
	MLPB93();

protected:
	/* Calculate the transition & transversion between two codons at a given position*/
	int TransitionTransversion(string codon1, string codon2, int pos);

}; 

#endif


