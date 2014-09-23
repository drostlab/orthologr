/*********************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
*
* Filename: LWL85.cpp
* Abstract: Definition of LWL85 and Modified LWL85 class.
*
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb., 2005

  References: 
  Li WH, Wu CI, Luo CC  (1985)  A new method for
  estimating synonymous and nonsynonymous rates of nucleotide 
  substitution considering the relative likelihood of nucleotide
  and codon changes. Mol. Biol. Evol. 2:150-174.

  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
**********************************************************/

#include "LWL85.h"


LWL85::LWL85() {
	int i;

	name = "LWL";	
	for(i=0; i<5; i++) K[i] = A[i] = B[i] = Pi[i] = Qi[i] = 0;
}

/************************************************
* Function: CountSiteAndDiff
* Input Parameter: codon1, codon2
* Output: Calculate synonymous and nonsynonymous
		 sites and differences between two codons.
* Return Value: void
*************************************************/
void LWL85::CountSiteAndDiff(string codon1, string codon2) {

	int i,j,k;
	int num=0, diff[5];
	string temp1 = "", temp2="";
	int sum=1, stop=0;

	//Count sites
	for(i=0; i<CODONLENGTH; i++) {		
		L[getCodonClass(codon1,i)]++;
		L[getCodonClass(codon2,i)]++;
	}
	
	//Count differences
	for (i=0; i<3; i++) {
		diff[i] = -1;
		if (codon1[i]!=codon2[i]) { 
			diff[num] = i;
			num++;
		}		
	}
	
	//Two codons are identical
	if(num==0) return;
	
	//Sum is the number of evolution pathway.
	for (i=1; i<=num; i++) sum *= i;

	snp += num;

	for(i=0; i<CODONLENGTH; i++)
		Si_temp[2*i] = Vi_temp[2*i] =0.0;
	
	//One difference
	if(num==1) {		
		TransitionTransversion(codon1, codon2, diff[0]);
	}

	//Two differences: 2 pathways for evolution (i,j)
	if(num==2) {
		for(i=0;i<num;i++) {
			for(j=0;j<num;j++) {
				if(i!=j) {
					temp1 = codon1;
					temp1[diff[i]] = codon2[diff[i]];
					if (getAminoAcid(temp1)!='!') {
						//pathway: codon1 <-> temp1 <-> codon2
						TransitionTransversion(codon1, temp1, diff[i]);
						TransitionTransversion(temp1, codon2, diff[j]);
					}
					else {
						stop++;
					}
				}
			}
		}	
	}

	//Three differences: 6 pathways for evolution (i,j,k)
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
							//pathway: codon1 <-> temp1 <-> temp2 <-> codon2
							TransitionTransversion(codon1, temp1, diff[i]);							
							TransitionTransversion(temp1,  temp2, diff[j]);														
							TransitionTransversion(temp2, codon2, diff[k]);
						}
						else
							stop++;
					}
				}
			}
		}
	}

	//Add pair-codon's differences to Si and Vi
	for(i=0; i<CODONLENGTH && (sum-stop)>0; i++) {
		Si[2*i] += (Si_temp[2*i]/(sum-stop));
		Vi[2*i] += (Vi_temp[2*i]/(sum-stop));
	}
}

/************************************************
* Function: getCodonClass
* Input Parameter: codon, position(0,1,2)
* Output: return 0,2,4-fold of codon at a given position.
* Return Value: 0 or 2 or 4
*************************************************/
int LWL85::getCodonClass(string codon, int pos) {
	int i;
	int codonClass = 0;

	for(i=0; i<4; i++) {
		if (i!=convertChar(codon[pos])) {
			string temp = codon;
			temp[pos] = convertInt(i);
			if (getAminoAcid(temp)!='!' && getAminoAcid(temp)==getAminoAcid(codon)) {
				codonClass++;
			}
		}
	}
	if (codonClass>0 && codonClass<3) {
		codonClass = 2;
	}
	else if (codonClass==3) {
		codonClass = 4;
	}

	return codonClass;		
}

/************************************************
* Function: TransitionTransversion
* Input Parameter: codon1, codon2, position(0,1,2)
* Output: Calculate synonymous and nonsynonymous differences
		of two codons at a given position.
* Return Value: int
* Note: Follow kakstools.py sent from Prof.Li
*************************************************/
int LWL85::TransitionTransversion(string codon1, string codon2, int pos) {
	
	//1:synonymous, 0:nonsynonymous
	int isSyn = 0;	

	/* Follow kakstools.py sent from Prof.Li WH */
	//CGA, CGG, AGA, AGG
	if ((codon1=="CGA" && codon2=="AGA" && pos==0)||(codon1=="AGA" && codon2=="CGA" && pos==0)) {
		isSyn = 1;
	}
	if ((codon1=="CGG" && codon2=="AGG" && pos==0)||(codon1=="AGG" && codon2=="CGG" && pos==0)) {
		isSyn = 1;
	}	
	if (isSyn==1) {		
		int c1=getCodonClass(codon1,pos);
		int c2=getCodonClass(codon2,pos);
		Si_temp[c1] += 0.5;	//different from the following 
		Vi_temp[c2] += 0.5;
				
		return 1;
	}

	/* Normal situation */
	//Synonymous: T(0)<->C(1), A(2)<->G(3) 	
	int sum = convertChar(codon1[pos]) + convertChar(codon2[pos]);
	if(sum==5 || sum==1)
		isSyn = 1;
	
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

	return 1;
}

/************************************************
* Function: preProcess
* Input Parameter: seq1, seq2
* Output: preprocess for Run
* Return Value: void
*************************************************/
void LWL85::preProcess(string seq1, string seq2) {

	long i;
	double ts=0, tv=0;
	double ai[5], bi[5];
	
	for(i=0; i<seq1.length(); i+=3) {		
		CountSiteAndDiff(seq1.substr(i,3), seq2.substr(i,3));
	}

	for(i=0; i<5; i+=2) {		
		
		ai[i] = bi[i] = 0.0;
		
		ts+=Si[i];
		tv+=Vi[i];

		L[i] = L[i]/2.0;
		Pi[i] = Si[i]/L[i];
		Qi[i] = Vi[i]/L[i];

		if ((1 - 2*Pi[i] - Qi[i])>0 &&  (1 - 2*Qi[i])>0) {
			ai[i] = 1/(1 - 2*Pi[i] - Qi[i]);
			bi[i] = 1/(1 - 2*Qi[i]);			
		}

		if (ai[i]>0 && bi[i]>0) {
			
			if (log(bi[i])>=0) {
				B[i] = 0.5*log(bi[i]);
			}
			
			if ((0.5*log(ai[i]) - 0.25*log(bi[i]))>=0) {
				A[i] = 0.5*log(ai[i]) - 0.25*log(bi[i]);
			}
			
			K[i] = A[i] + B[i];
		}

	}
	
	if(tv>0)
		kappa = 2*ts/tv;
	else
		kappa = 2;

	//For output formatting
	KAPPA[0] = KAPPA[1] = kappa;

	return;
}

/************************************************
* Function: Run
* Input Parameter: seq1, seq2
* Output: Main function for calculating Ka&Ks.
* Return Value: void
*************************************************/
string LWL85::Run(string seq1, string seq2) {
	
	preProcess(seq1, seq2);

	S = L[2]/3 + L[4];
	N = L[0] + 2*L[2]/3;
	
	Sd = L[2]*A[2] + L[4]*K[4];
	Ks = Sd / S;	
	
	Nd = L[0]*K[0] + L[2]*B[2];
	Ka = Nd/N;

	t = (S*Ks+N*Ka)/(S+N);

	return parseOutput();
}



MLWL85::MLWL85() {
	name = "MLWL";
}


/*For more detail see reference: Tzeng Y-H, Pan R, Li W-H  (2004)  Mol. Biol. Evol.*/
int MLWL85::TransitionTransversion(string codon1, string codon2, int pos) {
	
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

/* One of differences between MLWL85 and LWL85 is allowing for kappa in S and N */
string MLWL85::Run(string stra, string strb) {

	long i;
	double ts=0.0, tv=0.0;	//Transition, Transversion
	double ai[5], bi[5];	
	
	for(i=0; i<stra.length(); i+=3) {
		this->CountSiteAndDiff(stra.substr(i,3), strb.substr(i,3));
	}
	
	for(i=0; i<5; i+=2) {
		
		ai[i] = bi[i] = 0;

		ts+=Si[i];
		tv+=Vi[i];

		L[i] = L[i]/2.0;
		Pi[i] = Si[i]/L[i];
		Qi[i] = Vi[i]/L[i];

		if ((1 - 2*Pi[i] - Qi[i])>0 &&  (1 - 2*Qi[i])>0) {
			ai[i] = 1/(1 - 2*Pi[i] - Qi[i]);
			bi[i] = 1/(1 - 2*Qi[i]);			
		}
		
		if (ai[i]>0 && bi[i]>0) {
			
			if (log(bi[i])>=0) {
				B[i] = 0.5*log(bi[i]);
			}
			
			if ((0.5*log(ai[i]) - 0.25*log(bi[i]))>=0) {
				A[i] = 0.5*log(ai[i]) - 0.25*log(bi[i]);
			}
			
			K[i] = A[i] + B[i];
		}
	}
	
	kappa = 2.0*ts/tv;
	if (ts<SMALLVALUE || tv<SMALLVALUE) kappa = 1;

	KAPPA[0] = KAPPA[1] = kappa;

	if (kappa>2.0) {
		S = (kappa-1)*L[2]/(kappa+1) + L[4];
		N = L[0] + 2*L[2]/(kappa+1);		
	}
	else {
		if(kappa>0.5) {
			S = (kappa-0.5)*L[2]/(kappa+1.5) + L[4];
			N = L[0] + 2*L[2]/(kappa+1.5);		
		}
		else {
			S = L[2]/3 + L[4];
			N = 2*L[2]/3 + L[0];			
		}
	}

	Sd = L[2]*A[2] + L[4]*K[4];
	Ks = Sd/S;		

	Nd = L[0]*K[0] + L[2]*B[2];
	Ka = Nd/N;

	t = (S*Ks+N*Ka)/(S+N);
	
	return parseOutput();
}


