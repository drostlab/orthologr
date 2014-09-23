/************************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: base.cpp
* Abstract: Definition of base class for KaKs methods.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

*************************************************************/

#include "base.h"


/******** Global variables ********/


/*						The Genetic Codes 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
	Last update of the Genetic Codes: October 05, 2000 */
int genetic_code=1; //from 1 to 23
/* Genetic standard codon table, !=stop codon */
const char* transl_table[] = {
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "1-Standard Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS!!VVVVAAAADDEEGGGG", "2-Vertebrate Mitochondrial Code",
 "FFLLSSSSYY!!CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "3-Yeast Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "4-Mold Mitochondrial Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "5-Invertebrate Mitochondrial Code",
 "FFLLSSSSYYQQCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "6-Ciliate, Dasycladacean and Hexamita Code",
 "", "7-",
 "", "8-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "9-Echinoderm and Flatworm Mitochondrial Code",
 "FFLLSSSSYY!!CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "10-Euplotid Nuclear Code",
 "FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "11-Bacterial and Plant Plastid Code",
 "FFLLSSSSYY!!CC!WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "12-Alternative Yeast Nuclear Code",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "13-Ascidian Mitochondrial Code",
 "FFLLSSSSYYY!CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "14-Alternative Flatworm Mitochondrial Code",
 "FFLLSSSSYY!QCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "15-Blepharisma Nuclear Code",
 "FFLLSSSSYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "16-Chlorophycean Mitochondrial Code",
 "", "17-",
 "", "18-",
 "", "19-",
 "", "20-",
 "FFLLSSSSYY!!CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "21-Trematode Mitochondrial Code",
 "FFLLSS!SYY!LCC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "22-Scenedesmus obliquus mitochondrial Code",
 "FF!LSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "23-Thraustochytrium Mitochondrial Code"
};

string seq_name;		//sequences' name
unsigned long length;	//sequences' length
double GC[4];			//GC Content
//********End of Global variables**********

//Constructor function
Base::Base() {
	
	int i;
	for(i=0; i<5; i++) {
		Si[i] = Vi[i] = L[i] = NULL;
	}
	for (i=0; i<NUMBER_OF_RATES; i++) {
		KAPPA[i] = 1;
	}
	
	SEKa = SEKs = AICc = lnL = AkaikeWeight = NA;
	Ka = Ks = Sd = Nd = S = N = snp = t = kappa = NULL;

	model = "";
}

/********************************************
* Function: addString
* Input Parameter: string, string, string
* Output: result = result + str + flag
* Return Value: void
* Note: flag = "\t" (default) or "\n"
*********************************************/
void Base::addString(string &result, string str, string flag) {
	result += str;
	result += flag;
}

/**********************************************************************
* Function: getAminoAcid
* Input Parameter: codon or codon's id
* Output: Calculate the amino acid according to codon or codon's id.
* Return Value: char 
***********************************************************************/
char Base::getAminoAcid(string codon) {
	return transl_table[2*(genetic_code-1)][getID(codon)];
}
char Base::getAminoAcid(int id) {
	return transl_table[2*(genetic_code-1)][id];
}

/**********************************
* Function: getNumNonsense
* Input Parameter: int
* Output: get the number of nonsense codons
* Return Value: int
***********************************/
int Base::getNumNonsense(int genetic_code) {

	int num, i;
	for(num=i=0; i<CODON; i++) {
		if(getAminoAcid(i)=='!') num++;
	}

	return num;
}

/********************************************
* Function: getID
* Input Parameter: codon
* Output: Get codon's id in array of codon_table.
* Return Value: int
*********************************************/
int Base::getID(string codon) {
	return (convertChar(codon[0])*XSIZE + convertChar(codon[1])*DNASIZE + convertChar(codon[2]));
}

/********************************************
* Function: getCodon
* Input Parameter: int
* Output: Get the codon according to id;
		  a reverse funtion of getID.
* Return Value: string
*********************************************/
string Base::getCodon(int IDcodon) {
	
	string codon = "TTT";

	if (IDcodon>=0 && IDcodon<64) {
		codon[0]=convertInt(IDcodon/16); 
		codon[1]=convertInt((IDcodon%16)/4);
		codon[2]=convertInt(IDcodon%4);
	}

	return codon;
}

/*********************************************
* Function: convertChar
* Input Parameter: ch as char
* Output: Convert a char-T,C,A,G into a digit
*		  0,1,2,3, respectively.
* Return Value: int.
**********************************************/
int Base::convertChar(char ch) {
	int ret = -1;
	switch(ch) {
		case 'T':case 'U':
			ret = 0;
			break;
		case 'C':
			ret = 1;
			break;
		case 'A':
			ret = 2;
			break;
		case 'G':
			ret = 3;
			break;
	}
	return ret;
}

/********************************************
* Function: convertInt
* Input Parameter: int
* Output: Convert a digit- 0,1,2,3 into a 
*		  char-T,C,A,G, respectively.
* Return Value: char
*********************************************/
char Base::convertInt(int i) {
	char ch = '-';
	switch(i) {
		case 0:
			ch = 'T';
			break;
		case 1:
			ch = 'C';
			break;
		case 2:
			ch = 'A';
			break;
		case 3:
			ch = 'G';
			break;
	}
	return ch;
}

/********************************************
* Function: stringtoUpper
* Input Parameter: string
* Output: upper string
* Return Value: string
*********************************************/
string Base::stringtoUpper(string str) {
	int j;
	for(j=0; str[j] = toupper(str[j]), j<str.length(); j++);
	return str;
}

/********************************************
* Function: getRandom
* Input Parameter: void
* Output: Generate a radnom integer
* Return Value: int
*********************************************/
int Base::getRandom() {	
	srand((unsigned)time(NULL));
	return rand();
}

/********************************************
* Function: initArray
* Input Parameter: array of int/double, int, int/double(default=0)
* Output: Init the array x[0...n-1]=value
* Return Value: int
*********************************************/
int Base::initArray(double x[], int n, double value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

int Base::initArray(int x[], int n, int value) {
	int i; 
	for(i=0; i<n; i++) x[i] = value;
	return 0;
}

/********************************************
* Function: sumArray
* Input Parameter: double/int, int, int(default=0)
* Output: Sum of array x[]
* Return Value: double/int
*********************************************/
double Base::sumArray(double x[], int end, int begin) { 
	int i;
	double sum=0.;	
	for(i=begin; i<end; sum += x[i], i++);    
	return sum;
}

int Base::sumArray(int x[], int end, int begin) { 
	int i, sum=0;	
	for(i=begin; i<end; sum += x[i], i++);    
	return sum;
}

/********************************************
* Function: norm
* Input Parameter: array of double, int
* Output: Sqrt of the sum of the elements' square 
           sqrt(x0*x0 + x1*x1 + ...)
* Return Value: double
*********************************************/
double Base::norm(double x[], int n) {
	int i; 
	double t=0; 

	for(i=0; i<n; t += square(x[i]), i++);

	return sqrt(t);
}

/********************************************
* Function: scaleArray
* Input Parameter: double, array of double, int
* Output: Elements in array are mutipled by scale 
* Return Value: int
*********************************************/
int Base::scaleArray(double scale, double x[], int n) {	
	int i; 	
	for (i=0; i<n; i++) x[i] *= scale;

	return 1; 
}

/********************************************
* Function: copyArray
* Input Parameter: array, array, int
* Output: Copy array's values one by one: to[] = from[]
* Return Value: int
*********************************************/
int Base::copyArray(double from[], double to[], int n) {	
	int i; 
	for(i=0; i<n; i++) to[i] = from[i];
	
	return 1; 
}

/********************************************
* Function: innerp
* Input Parameter: array, array, int
* Output: Sum of 'n' products multiplied by 
			two elements in x[], y[].
* Return Value: int
*********************************************/
double Base::innerp(double x[], double y[], int n) {
	
	int i; 
	double t=0;

	for(i=0; i<n; t += x[i]*y[i], i++); 

	return t; 
}

/********************************************
* Function: initIdentityMatrix
* Input Parameter: array of double, int
* Output: Set x[i,j]=0 when x!=j and 
			  x[i,j]=1 when x=j 
* Return Value: int
*********************************************/
int Base::initIdentityMatrix(double x[], int n) {
	
	int i,j;

	for (i=0; i<n; i++)  {
		for(j=0; j<n; x[i*n+j]=0, j++);
		x[i*n+i] = 1; 
	} 

	return 0; 
}


/************************************************
* Function: writeFile
* Input Parameter: string, string
* Output: Write content into the given file.
* Return Value: True if succeed, otherwise false. 
*************************************************/
bool Base::writeFile(string output_filename, const char* result) {
	
	bool flag = true;
	try {
		//file name is ok
		if (output_filename!="" && output_filename.length()>0) {
			ofstream os(output_filename.c_str());
			if (!os.is_open()) throw 1;

			os<<result;
			os.close();					
		}
	}
	catch (...) {
		cout<<"Error in writing to file..."<<endl;
		flag = false;
	}	

	return flag;
}


/*****************************************************
* Function: parseOutput
* Input Parameter: void
* Output: Parse estimated results for outputing
* Return Value: string

  Order: "Sequence", "Method", "Ka", "Ks", "Ka/Ks", 
		 "P-Value(Fisher)", "Length", "S-Sites", "N-Sites", "Fold-Sites(0:2:4)",
		 "Substitutions", "S-Substitutions", "N-Substitutions", "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", 
		 "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
		 "Model"
******************************************************/		
string Base::parseOutput() {

	int i;
	string result = "", tmp;

	//Sequence name
	addString(result, seq_name);
	//Method name
	addString(result, name);

	//Ka
	if (Ka<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ka);
	}
	addString(result, tmp);
	
	//Ks
	if (Ks<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ks);
	}
	addString(result, tmp);
	
	//Ka/Ks
	if(Ks<SMALLVALUE || Ks==NA || Ka==NA) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(Ka/Ks);
	}
	addString(result, tmp);

	//Fisher's test: p_value
	if (Sd<SMALLVALUE || Nd<SMALLVALUE || S<SMALLVALUE || N<SMALLVALUE)
		tmp = "NA";
	else {
		tmp = CONVERT<string>(fisher(Sd,Nd,S-Sd,N-Nd));
	}
	addString(result, tmp);

	//Length of compared pairwise sequences
	addString(result, CONVERT<string>(length));
	
	//Synonymous(S) sites
	if (S<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(S);
	}
	addString(result, tmp);

	//Nonsynonymous(N) sites
	if (N<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(N);
	}
	addString(result, tmp);

	//L[0], L[2], L[4] only for Prof.Li's series(LWL85, LPB93...)
	if (L[0]<SMALLVALUE && L[2]<SMALLVALUE && L[4]<SMALLVALUE) {
		tmp = "NA";		
	}
	else {		
		tmp = CONVERT<string>(L[0]);	tmp += ":";
		tmp += CONVERT<string>(L[2]);	tmp += ":";
		tmp += CONVERT<string>(L[4]);		
	}
	addString(result, tmp);

	
	//Substitutions
	addString(result, CONVERT<string>(snp));

	//Sysnonymous(Sd) Substitutions(Nd)	
	if (Sd>SMALLVALUE) {
		tmp = CONVERT<string>(Sd);		
	}	
	else {
		tmp = "NA";
	}
	addString(result, tmp);
	
	//Nonsysnonymous Substitutions(Nd)
	if (Nd>SMALLVALUE) {
		tmp = CONVERT<string>(Nd);		
	}	
	else {
		tmp = "NA";
	}
	addString(result, tmp);

	//Si for Li's series' methods(LWL85, LPB93...)
	if (Si[0]!=NULL || Si[2]!=NULL || Si[4]!=NULL) { //Si[0], Si[2], Si[4]
		tmp  = CONVERT<string>(Si[0]);	tmp += ":";
		tmp += CONVERT<string>(Si[2]);	tmp += ":";
		tmp += CONVERT<string>(Si[4]);		
	}
	else {
		tmp = "NA";
	}
	addString(result, tmp);

	//Vi for Li's series' methods(LWL85, LPB93...)
	if (Vi[0]!=NULL || Vi[2]!=NULL || Vi[4]!=NULL) { //Vi[0], Vi[2], Vi[4]
		tmp  = CONVERT<string>(Vi[0]);	tmp += ":";
		tmp += CONVERT<string>(Vi[2]);	tmp += ":";
		tmp += CONVERT<string>(Vi[4]);		
	}
	else {
		tmp = "NA";
	}
	addString(result, tmp);


	//Divergence time or distance t = (S*Ks+N*Ka)/(S+N)
	if(t<SMALLVALUE) {
		tmp = "NA";
	}
	else {
		tmp = CONVERT<string>(t);
	}
	addString(result, tmp);

	//Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)
	for(i=0, tmp=""; i<NUMBER_OF_RATES-1; i++) {
		tmp += CONVERT<string>(KAPPA[i]); 
		tmp += ":";
	}
	tmp += CONVERT<string>(KAPPA[i]); 
	addString(result, tmp);

	//GC Content
	tmp  = CONVERT<string>(GC[0]);	tmp += "(";
	tmp += CONVERT<string>(GC[1]);	tmp += ":";
	tmp += CONVERT<string>(GC[2]);	tmp += ":";
	tmp += CONVERT<string>(GC[3]);	tmp += ")";
	addString(result, tmp);
	
	//Maximum Likelihood Value
	if (lnL==NA) tmp = "NA";
	else         tmp = CONVERT<string>(lnL);
	addString(result, tmp);
	
	//AICc
	if (AICc==NA) tmp = "NA";
	else		  tmp = CONVERT<string>(AICc);
	addString(result, tmp);

	//Akaike weight in model selection
	if (AkaikeWeight==NA) tmp = "NA";
	else		  tmp = CONVERT<string>(AkaikeWeight);
	addString(result, tmp);
	
	//Selected Model according to AICc
	if (model==""||model.length()==0) tmp = "NA";
	else tmp = model;
	addString(result, tmp, "\n");

/*
	//Standard Errors
	if (SEKa==NA) tmp = "NA";
	else tmp = CONVERT<string>(SEKa);
	addString(result, tmp, "\t");
	
	if (SEKs==NA) tmp = "NA";
	else tmp = CONVERT<string>(SEKs);
	addString(result, tmp, "\n");
*/
	
	return result;
}

/**************************************************
* Function: fisher
* Input Parameter: double, double, double, double
* Output: Compute p-value by Fisher exact test
* Return Value: double
***************************************************/
double Base::fisher(double sd, double nd, double s, double n) {
	double denominator, numerator, prob_total, prob_current, sum, fac_sum;
	double matrix[4],  R[2], C[2];
	int i, j;

	denominator = numerator = prob_total = prob_current = sum = fac_sum = 0.0;

	matrix[0]=sd;      matrix[1]=s;	    matrix[2]=nd;    matrix[3]=n;
	//Row & Column
	R[0]=matrix[0]+matrix[2];		R[1]=matrix[1]+matrix[3];
	C[0]=matrix[0]+matrix[1];		C[1]=matrix[2]+matrix[3];
	sum = R[0]+R[1];	
	
	//Calculate the numberator that is a constant
	numerator += factorial(R[0]);
	numerator += factorial(R[1]);
	numerator += factorial(C[0]);
	numerator += factorial(C[1]);
	
	//Log of Factorial of N
	fac_sum = factorial(sum);	
	for(i=0, denominator=fac_sum; i<4; i++) {
		denominator += factorial(matrix[i]);
	}
	//Probability of current situtation
	prob_current = exp(numerator-denominator);
	
	//Two-tail probabilities if less than prob_current
	for(i=0; i<=R[0]; i++) {
		matrix[0] = i;       
		matrix[1] = C[0]-i;
		matrix[2] = R[0]-i;
		matrix[3] = R[1]-C[0]+i;
		if (matrix[0]>=0 && matrix[1]>=0 && matrix[2]>=0 && matrix[3]>=0) {			
			for(j=0, denominator=fac_sum; j<4; j++) {				
				denominator += factorial(matrix[j]);
			}
			double temp = numerator-denominator;
			temp = exp(numerator-denominator);
			if (temp<=prob_current) {
				prob_total += temp;
			}
		}
	}
	
	return prob_total;
}

/**************************************************
* Function: factorial
* Input Parameter: n
* Output: Compute the factorial of 'n', then return 
          the log of it.
* Return Value: double
***************************************************/
double Base::factorial(double n) {
	double temp=1.0;
	if (n>0) {
		n = n + 1;
		double x = 0;
		x += 0.1659470187408462e-06/(n+7);
		x += 0.9934937113930748e-05/(n+6);
		x -= 0.1385710331296526    /(n+5);
		x += 12.50734324009056     /(n+4);
		x -= 176.6150291498386     /(n+3);
		x += 771.3234287757674     /(n+2);
		x -= 1259.139216722289     /(n+1);
		x += 676.5203681218835     /(n);
		x += 0.9999999999995183;
		temp = log(x) - 5.58106146679532777 - n +(n - 0.5)*log(n + 6.5);
	}
	
	return (temp);
}


/*

int matby (double a[], double b[], double c[], int n,int m,int k)
// a[n*m], b[m*k], c[n*k]  ......  c = a*b

{
   int i,j,i1;
   double t;
   FOR (i,n)  FOR(j,k) {
      for (i1=0,t=0; i1<m; i1++) t+=a[i*m+i1]*b[i1*k+j];
      c[i*k+j] = t;
   }
   return (0);
}





void PMatrixTaylor(double P[], double n, double t) {

// Get approximate PMat using polynomial approximation
//   P(t) = I + Qt + (Qt)^2/2 + (Qt)^3/3!

   int nterms=1000, i,j,k, c[2],ndiff,pos=0,from[3],to[3];
   double *Q=space, *Q1=Q+n*n, *Q2=Q1+n*n, mr, div;

   FOR (i,n*n) Q[i]=0;
   for (i=0; i<n; i++) FOR (j,i) {
      from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
      to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
      c[0]=GenetCode[com.icode][i];
      c[1]=GenetCode[com.icode][j];
      if (c[0]==-1 || c[1]==-1)  continue;
      for (k=0,ndiff=0; k<3; k++)  if (from[k]!=to[k]) { ndiff++; pos=k; }
      if (ndiff!=1)  continue;
      Q[i*n+j]=1;
      if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  Q[i*n+j]*=kappa;
      if(c[0]!=c[1])  Q[i*n+j]*=omega;
      Q[j*n+i]=Q[i*n+j];
   }
   FOR(i,n) FOR(j,n) Q[i*n+j]*=pi[j];
   for (i=0,mr=0;i<n;i++) { 
      Q[i*n+i]=-sum(Q+i*n,n); mr-=pi[i]*Q[i*n+i]; 
   }
   FOR(i,n*n) Q[i]*=t/mr;

   xtoy(Q,P,n*n);  FOR(i,n) P[i*n+i]++;   // I + Qt 
   xtoy(Q,Q1,n*n);
   for (i=2,div=2;i<nterms;i++,div*=i) {  // k is divisor 
      matby(Q, Q1, Q2, n, n, n);
      for(j=0,mr=0;j<n*n;j++) { P[j]+=Q2[j]/div; mr+=fabs(Q2[j]); }
      mr/=div;
      // if(debug) printf("Pmat term %d: norm=%.6e\n", i, mr); 
      if (mr<e) break;
      xtoy(Q2,Q1,n*n);
   }

   FOR(i,n*n) if(P[i]<0) P[i]=0;

}

*/
