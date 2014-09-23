/************************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: GY94.h
* Abstract: Declaration of GY94 class.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Oct.2, 2005

* Note: Source codes are taken from codeml.c in PAML.

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.
*************************************************************/

#if !defined(GY94_H)
#define  GY94_H

#include "base.h"

using namespace std;

class GY94: public Base {

public:	
	GY94();
	GY94(string model);
	~GY94();
	
	/* Main function */
	string Run(const char *seq1, const char *seq2);
	
protected:
	/* Preprocess for calculating Ka&Ks */
	int preProcess(const char* seq1, const char* seq2);
	/* Parse substitution rates according to the given model */
	int parseSubRates(string model, double kappa[]);
	/* Construct two array according to genetic code */
	int setmark_61_64 (void);	
	int PatternWeight();
	/* Encode two compared sequneces */
	void EncodeSeqs (void);	
	/* Calculate Ka&Ks using ML */
	int PairwiseCodon (double space[]);
	/* Get codons' frequencies */
	int GetCodonFreqs(double pi[]);	
	/* Construct transition probability matrix (64*64) */
	int EigenQc (int getstats, double blength, double *S, double *dS, double *dN, double Root[], double U[], double V[], double kappa[], double omega, double Q[]);
	/* Return maximum-likelihood score */
	double lfun2dSdN (double x[], int np);
	/* Main fuctiion for GY method */
	int ming2 (double *f, double x[], double xb[][2], double space[], double e, int n);	

	int eigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[]);	
	double LineSearch2 (double *f, double x0[], double p[], double step, double limit, double e, double space[], int n);
	double rndu(void) ;
	/* x[i]=x0[i] + t*p[i] */
	double fun_ls(double t, double x0[], double p[], double x[], int n);	
	double distance (double x[], double y[], int n);		
	int transform (char *z, int ls);	
	int gradientB (int n, double x[], double f0, double g[], double space[], int xmark[]);	
	int H_end (double x0[], double x1[], double f0, double f1, double e1, double e2, int n);
	void HouseholderRealSym(double a[], int n, double d[], double e[]);
	int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
	void EigenSort(double d[], double U[], int n);
	int eigenRealSym(double A[], int n, double Root[], double work[]);
	void FreeMemPUVR(void);

public:
	struct common_info {
		char *z[2];	//sequences
		int ns;	//ns=2
		int ls;	//sequence's length
		int ngene, npatt; 
		int icode;	//number of genetic code
		int ncode;	//number of non-stop codon
		int np, nkappa, sspace;
		double fpatt[CODON*CODON]; 
		double *space;
		double kappa;	//transition/transversion
		double omega;	//Ka/Ks
		double pi[CODON];	//codons' frequencies
		double KAPPA[8]; 
	}  com;

protected:
	double PMat[CODON*CODON],U[CODON*CODON],V[CODON*CODON],Root[CODON*CODON];
	int Nsensecodon, FROM61[CODON], FROM64[CODON];
	int Iround;
	unsigned int w_rndu;//=123456757;
	double SIZEp;
	double Small_Diff; 
	string str1, str2;	//a pair of sequences
};

#endif
