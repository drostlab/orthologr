/******************************************************************
* Copyright (c) 2006, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: base.h
* Abstract: Declaration of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

* Modified Version: 1.2 
* Modified Author: ZZ
* Modified Date: Apr. 2006
******************************************************************/
#pragma warning(disable:4786)		//Disable warning when using vector

#if !defined(BASE_H)
#define  BASE_H

#define VERSION		"1.2, Apr. 2006"

#define CODONLENGTH 3			//Length of codon
#define DNASIZE 4				//A C G T
#define XSIZE DNASIZE*DNASIZE   //Size of one group AXX (X=A,C,G,T) 
#define CODON 64				//Codon Size
#define NULL 0					//Zero
#define NA -1					//Not Available
#define NCODE	23				//Number of genetic codes
#define NNCODE  NCODE*2			//Double of the number genetic codes		
#define SMALLVALUE 1e-6			//Value near to zero
#define NUMBER_OF_RATES	6		//Number of substitution rates
#define MODELCOUNT	14			//Number of candidate models

#define gammap(x,alpha) (alpha*(1-pow(x,-1/alpha)))
#define square(a) ((a)*(a))

#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Stanard lib of C++ */
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<time.h>


using namespace std;

/*****Global variables*****/
extern const char* transl_table[];	//Genetic codon table
extern int genetic_code;			//ID of codon table from 1 to 23

extern string seq_name;				//Pairwise sequences' name
extern unsigned long length;		//Length of compared sequences
extern double GC[4];				//GC Content of entire sequences(GC[0]) and three codon positions (GC[1--3])
//End of Global variables


/* Convert one type to any other type */
template<class out_type,class in_value>
	out_type CONVERT(const in_value & t) {
		stringstream stream;
		//Put the value 't' into the stream
		stream<<t;			
		out_type result;
		//Put the stream into the 'result'
		stream>>result;

		return result;
	}

class Base {

public:
	Base();

	/* Write the content into file */
	bool writeFile(string output_filename, const char* result);

	/* Parse results */
	string parseOutput();	
	/* Format string for outputing into file */
	void addString(string &result, string str, string flag="\t");

	/* Generate a radnom integer */
	int getRandom();
	
	/* Convert a char-T,C,A,G into a digit 0,1,2,3, respectively */
	int  convertChar(char ch);
	/* Convert a digit-0,1,2,3 into a char T,C,A,G, respectively */
	char convertInt(int ch);
	/* Convert a string to uppercase */
	string stringtoUpper(string str);
	
	/* Calculate the amino acid of codon by codon */
	char getAminoAcid(string codon);
	/* Calculate the amino acid of codon by codon's id*/
	char getAminoAcid(int id);
	/* Get the number of stop codon in a given genetic code table */
	int getNumNonsense(int genetic_code);

	/* Return the codon's id from codon table */
	int getID(string codon);
	/* Return a codon according to the id */
	string getCodon(int IDcodon);

	/* Sum array's elements */
	double sumArray(double x[], int end, int begin=0);
	int sumArray(int x[], int end, int begin=0);

	/* Init value to array */
	int initArray(double x[], int n, double value=0.0);
	int initArray(int x[], int n, int value=0);

	/* Elements in array are mutipled by scale */
	int scaleArray(double scale, double x[], int n);
	/* Sqrt of the sum of the elements' square */
	double norm(double x[], int n);
	/* Copy array's values one by one: to[] = from[] */
	int copyArray(double from[], double to[], int n);
	/* Sum of 'n' products multiplied by two elements x[], y[] */
	double innerp(double x[], double y[], int n);
	/* Set x[i,j]=0 when x!=j and x[i,j]=1 when x=j */
	int initIdentityMatrix(double x[], int n);

protected:
	/* Compute p-value by Fisher exact test to justify the validity of ka/ks */
	double fisher(double cs, double us, double cn, double un);
	/* factorial */
	double factorial(double n);


protected:
	/* Name of method for calculating ka/ks */
	string name;
	/* Sysnonymous sites(S) and nonsysnonymous sites(N) */
	double S, N;
	/* Number of sysnonymous differences(Sd) and nonsysnonymous differences(Nd), snp=Sd+Nd */
	double Sd, Nd, snp;
	/* Number of sysnonymous substitutions(ks) and nonsysnonymous substitutions(ka) per site */
	double Ka, Ks;
	/* Standard Errors for Ka and Ks */
	double SEKa, SEKs;
	/* Transitional(Si) and transversional(Vi) substitutions */
	double Si[5], Vi[5];
	/* 0-fold, 2-fold, 4-fold */
	double L[5];
	/* Total Numbers of substitutions per i-th type site: K[i]=A[i]+B[i] */
	double K[5];
	/* 	   T  C  A  G
		T  -  6  7  8  
		C  0  -  9  10
		A  1  3  -  11
		G  2  4  5  -				*/	

	/* Transition/transversion rate ratio, tc: between pyrimidines, ag: between purines */
	double kappa, kappatc, kappaag;

public:
	/* Divergence time, substitutions per site: t = (S*Ks+N*Ka)/(S+N)*/
	double t;
	/* Maximum Likelihood value */
	double lnL;
	/* Akaike Information Criterion (AICc) */
	double AICc;
	/* Akaike Weight in Model Selection */
	double AkaikeWeight;
	/* Name of mutation model according to AICc */
	string model;
	/* Transition/transversion mutation ratio(s) or substitution rates */
	double KAPPA[NUMBER_OF_RATES];

protected:
	/* Store Maximum Likilhood results for Model Selection and Model Averaging */
	struct MLResult {
		string result;	//estimated results
		double AICc;	//the value of a modified AIC
		double freq[CODON];		//Codon frequency
		double rate[NUMBER_OF_RATES];	//Six subsitution rates
		double w;	//Ka/Ks
		double t;	//divergence time
	};
	
};

#endif


