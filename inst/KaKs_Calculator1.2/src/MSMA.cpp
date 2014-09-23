/************************************************************
* Copyright (c) 2006, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: MSMA.cpp
* Abstract: Definition of model-selected and model-averged methods' classes.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Apr. 2006

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.

  Posada, D. and Buckley, T.R. (2004) Model Selection and Model Averaging
  in Phylogenetics: Advantages of Akaike Information Criterion and Bayesian
  Approaches over Likelihood Ratio Tests, Syst. Biol., 53, 793-808.

  Sullivan, J. and Joyce, P. (2005) Model Selection in Phylogenetics, 
  Annual Review of Ecology, Evolution, and Systematics, 36, 445-466.  
*************************************************************/
#include "MSMA.h"

MS::MS() {
	name = "MS";
}

/* Calculate Ka and Ks based on a given model, similar to the method of GY */
void MS::selectModel(const char *seq1, const char *seq2, string c_model, vector<MLResult>& result4MA) {

	MLResult tmp;
	GY94 zz(c_model);	

	tmp.result = zz.Run(seq1, seq2);
	tmp.AICc = zz.AICc;
	copyArray(zz.com.pi, tmp.freq, (int)CODON);
	copyArray(zz.KAPPA, tmp.rate, (int)NUMBER_OF_RATES);
	tmp.w = zz.com.omega;
	tmp.t = 3.*zz.t;

	result4MA.push_back(tmp);
}

/* Choose the estimates under a model with smallest AICc */
string MS::Run(const char *seq1, const char *seq2, vector<MLResult> &result4MA, string &details) {

	int i, j, pos;
	string candidate_models[] = {"JC", "F81", "K2P", "HKY", "TNEF", "TN", "K3P", "K3PUF", "TIMEF", "TIM", "TVMEF", "TVM", "SYM", "GTR"};
	
	//Calculate Ka and Ks using 14 models
	for (i=0; i<MODELCOUNT; i++) selectModel(seq1, seq2, candidate_models[i], result4MA);

	//Choose the results under a model with smallest AICc
	for (pos=i=0; i<result4MA.size(); i++) {		
		if (result4MA[i].AICc<result4MA[pos].AICc) pos = i;
	}
	
	//Calculate the AICc difference, substract the smallest AICc
	double diff[MODELCOUNT];
	for (i=0; i<MODELCOUNT; i++) diff[i] = result4MA[i].AICc - result4MA[pos].AICc;

	//Compute Akaike weights
	double w[MODELCOUNT], sum;
	initArray(w, MODELCOUNT);
	for (sum=i=0; i<MODELCOUNT; i++) {
		//akaike wights of each model
		for (j=0; j<MODELCOUNT; j++) {
			double power = -0.5*diff[j] - (-0.5*diff[i]);
			//Avoid overflow
			if (power>709) power = 700;
			else if (power<-709) power = -700;
			//Normal
			w[i] += exp(power);
		}
		w[i] = 1./w[i];
	}

	//Add Akaike weights to results
	string tmp="";
	for (i=0; i<MODELCOUNT; i++) {

		//Add akaike weights
		tmp = result4MA[i].result;
		j = tmp.find_last_of('\t');
		result4MA[i].result = tmp.substr(j, tmp.length()-j);
		tmp = tmp.replace(j, tmp.length()-j, "");
		j = tmp.find_last_of('\t');
		tmp = tmp.replace(j+1, tmp.length()-j-1, "") + CONVERT<string>(w[i]);
		result4MA[i].result = tmp + result4MA[i].result;

		//Details on model selection
		details += result4MA[i].result;
	}

	//Results at "pos" is more reliable, replace 'method name' by "MS".
	tmp = result4MA[pos].result;	
	i = tmp.find('\t');
	j = tmp.find('\t', i+1);	
	result4MA[pos].result  = tmp.substr(0, i+1) + name;
	result4MA[pos].result += tmp.substr(j, tmp.length()-j);

	return result4MA[pos].result;
}


MA::MA() {
	name = "MA";
	Small_Diff=1e-6; 
	w_rndu=123456757;
	com.ns = 2; 
	lnL = 0.0;

	com.icode = genetic_code-1;
	if (com.icode>11) com.icode = 0;	
	com.ncode = Nsensecodon = 64 - getNumNonsense(com.icode);

	com.nkappa = 5;
	com.np = 2 + com.nkappa;
}


string MA::Run(const char *seq1, const char *seq2, vector<MLResult> result4MA) {

	int i, j, pos;
	
	//Find the smallest AICc
	for (pos=0, i=1; i<result4MA.size(); i++) {
		if (result4MA[i].AICc<result4MA[pos].AICc) pos = i;
	}

	//Calculate the AICc difference, substract the smallest AICc
	double diff[MODELCOUNT];
	for (i=0; i<MODELCOUNT; i++) diff[i] = result4MA[i].AICc - result4MA[pos].AICc;

	//Compute Akake weights
	double w[MODELCOUNT];
	initArray(w, MODELCOUNT);
	for (i=0; i<MODELCOUNT; i++) {
		//Avoid overflow
		for (j=0; j<MODELCOUNT; j++) {
			double power = -0.5*diff[j] - (-0.5*diff[i]);
			if (power>709) power = 700;
			else if (power<-709) power = -700;
			w[i] += exp(power);
		}
		w[i] = 1./w[i];
	}
	
	int I[MODELCOUNT][NUMBER_OF_RATES];
	for(i=0; i<MODELCOUNT; i++) 
		for (j=0; j<NUMBER_OF_RATES; j++)
			I[i][j]=0;
	
	//JC, F81:	  rTC==rAG =rTA==rCG==rTG==rCA

	//K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
	I[2][0]=1;	I[2][1]=1;  
	I[3][0]=1;	I[3][1]=1;
	//TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
	I[4][0]=1;	I[4][1]=1;
	I[5][0]=1;	I[5][1]=1;
	//K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
	I[6][0]=1;	I[6][1]=1;  I[6][2]=1;	I[6][3]=1;
	I[7][0]=1;	I[7][1]=1;	I[7][2]=1;	I[7][3]=1;
	//TIMEF, TIM: rTC!=rAG!=rTA==rCG!=rTG==rCA
	I[8][0]=1;	I[8][1]=1;  I[8][2]=1;	I[8][3]=1;
	I[9][0]=1;	I[9][1]=1;	I[9][2]=1;	I[9][3]=1;
	//TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
	I[10][0]=1;	I[10][1]=1;  I[10][2]=1;	I[10][3]=1;	I[10][4]=1;
	I[11][0]=1;	I[11][1]=1;	 I[11][2]=1;	I[11][3]=1;	I[11][4]=1;
	//SYM, GTR:   rTC!=rAG!=rTA!=rCG!=rTG!=rCA
	I[12][0]=1;	I[12][1]=1;  I[12][2]=1;	I[12][3]=1;	I[12][4]=1;
	I[13][0]=1;	I[13][1]=1;	 I[13][2]=1;	I[13][3]=1;	I[13][4]=1;

		
	//Average parameters
	double para[8];
	initArray(para, 8);
	
	//Model-averaged time t
	for(j=0; j<MODELCOUNT;j++) para[0] += (w[j]*result4MA[j].t);
	//para[0] /= sum;

	//Model-averaged Substitution Rates
	double sum[NUMBER_OF_RATES];
	initArray(sum, NUMBER_OF_RATES);
	initArray(com.KAPPA, 8);

	for (i=0; i<NUMBER_OF_RATES-1; i++) {
		for(j=0; j<MODELCOUNT;j++) {			
			com.KAPPA[i] += (w[j]*I[j][i]*result4MA[j].rate[i]);
			sum[i] += (w[j]*I[j][i]);			
		}
		if(sum[i]<1e-50) com.KAPPA[i]=1;
		else com.KAPPA[i] /= sum[i];
		para[i+1] = com.KAPPA[i];
	}

	para[6] = 0;
	//Model-averaged omega w
	for(j=0; j<MODELCOUNT;j++) para[6] += (w[j]*result4MA[j].w);
	//para[6] /= sum;

	//Model-averaged Codon Frequencies
	for (i=0; i<CODON; i++) {
		com.pi[i] = 0.0;
		for(j=0; j<MODELCOUNT;j++) com.pi[i] += (w[j]*result4MA[j].freq[i]);
		//com.pi[i] /= sum;
	}

	/* Preprocess */
	preProcess(seq1, seq2);
	
	//Calculate maximum likelihood score
	lnL = lfun2dSdN(para, com.np);

	//Copy subsitution rates
	copyArray(com.KAPPA, KAPPA, NUMBER_OF_RATES);
	
	//Compute Ka and Ks
	EigenQc(1, para[0], &S, &Ks, &Ka, NULL, NULL, NULL, KAPPA, com.omega, PMat);

	N = com.ls*3-S;
	Sd = Ks*S;
	Nd = Ka*N;	
	double scale=(Sd+Nd)/snp;
	Sd /= scale;
	Nd /= scale;

	lnL = -lnL;

	t = para[0]/3;
	
	//For the method of Model averaging, parameters' number is not specific.
	AICc = NA;	

	return parseOutput();

}