/**********************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: YN00.h
* Abstract: Declaration of method YN00.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

* Note: Source codes are taken from yn00.c in PAML.

  References:
  Yang Z, Nielsen R  (2000)  Estimating Synonymous and 
  Nonsynonymous Substitution Rates Under Realistic 
  Evolutionary Models. Mol Biol Evol 17:32-43.
**********************************************************/

#if !defined(YN00_H)
#define  YN00_H

#define CODONFREQ 12
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define square(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#include "base.h"

using namespace std;

class YN00: public Base {

public:
	YN00();	

	/* Main function of calculating kaks */
	string Run(string seq1, string seq2);
	
//protected:
	/* Get A,C,G,T's frequency between pair sequences: f12pos[], pi[], pi_sqrt[]  */
	void getFreqency(const string seq1, const string seq2);
	/* Get the k(transition/transversion) */
	virtual int GetKappa(const string seq1, const string seq2);
	/* Use the HKY85 Model to correct for multiple substitutions */
	virtual int DistanceF84(double n, double P, double Q, double pi4[],double &k_HKY, double &t, double &SEt);
	/* Calculate the ka,ks */
	virtual int DistanceYN00(const string seq1, const string seq2, double &dS,double &dN, double &SEdS, double &SEdN);
	/* Count synonymous and nonsynonmous sites: S, N */
	virtual int CountSites(const string z, double &Stot,double &Ntot,double fbS[],double fbN[]);
	/* Calculate the transition probability matrix using 'kappa' and 'omega' */
	virtual int GetPMatCodon(double P[], double kappa, double omega);
	/* Count synonymous and nonsynonmous differences: Sd, Nd */	
	virtual int CountDiffs(const string seq1, const string seq2, double &Sdts,double &Sdtv,double &Ndts, double &Ndtv,double PMatrix[]);
	
	//The following is for calculation of transition probability matrix by Taylor equation
	int eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0, double Root[], double U[], double V[]);
	int PMatUVRoot (double P[], double t, int n, double U[], double V[], double Root[]);
	int eigenRealSym(double A[], int n, double Root[], double work[]);
	void EigenSort(double d[], double U[], int n);
	int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
	void HouseholderRealSym(double a[], int n, double d[], double e[]);
	
public:
	double omega;	//Ka/Ks

protected:	
	/* 	  T	C A G
		1 *	* * *  
		2 *	* * *
		3 *	* * *
	f12pos: 3*4 matrix*/
	/* Probability of A,C,G,T at three positions */
	double f12pos[CODONFREQ];		
	/* Probability of 64 codons */
	double pi[CODON];
	double pi_sqrt[CODON];
	/* Whether iteration calculating ka/ks */
	int iteration;
	/* The number of which pi[] is zero */
	int npi0;

};

#endif

