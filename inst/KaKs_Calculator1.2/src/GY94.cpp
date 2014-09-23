/************************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: GY94.cpp
* Abstract: Definition of GY94 class.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Oct.2, 2005

* Note: Source codes are taken from codeml.c in PAML.

  References:
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.
*************************************************************/
#include "GY94.h"
//#include <malloc.h>
//#include<MALLOC>

int GeneticCode[][64] = 
     {{13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 0:universal */

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15,-1,-1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 1:vertebrate mt.*/

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       16,16,16,16,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 2:yeast mt. */

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 3:mold mt. */

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15,15,15,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 4:invertebrate mt. */

      {13,13,10,10,15,15,15,15,18,18, 5, 5, 4, 4,-1,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 5:ciliate nuclear*/

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2, 2,11,15,15,15,15,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 6:echinoderm mt.*/

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4, 4,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 7:euplotid mt. */

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
       10,10,10,15,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},
                                                 /* 8:alternative yeast nu.*/

      {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15, 7, 7,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 9:ascidian mt. */

      {13,13,10,10,15,15,15,15,18,18,-1, 5, 4, 4,-1,17,
       10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
        9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
       19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}, /* 10:blepharisma nu.*/

      { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
        5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
        9, 9, 9, 9,10,10,10,10,11,11,11,11,12,12,12,12,
       13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16} /* 11:Ziheng's regular code */
     };                                         /* GeneticCode[icode][#codon] */


/********************************************
* Function: construction function
* Input Parameter: a candidate model
* Output: 
* Return Value:

* Note: (Equal, Unequal: substitution rates)
JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
*********************************************/

//Constructor
GY94::GY94() {
}

GY94::GY94(string NulModel) {

	name = "GY-"+NulModel;
	Iround=0;
	SIZEp=0;
	Small_Diff=1e-6; 
	w_rndu=123456757;
	com.ns = 2; 
	lnL = 0.0;
	
	com.icode = genetic_code-1;
	if (com.icode>11) com.icode = 0;
	Nsensecodon = getNumNonsense(com.icode);
	com.ncode = 64 - Nsensecodon;
	
	//NulModel is set among a set of candidate models.
	model = NulModel;
	
	if (model=="JC" || model=="F81") com.nkappa = 0;
	else if (model=="K2P" || model=="HKY") com.nkappa = 1;
	else if (model=="TNEF"|| model=="TN") com.nkappa = 2;
	else if (model=="K3P" || model=="K3PUF") com.nkappa = 2;
	else if (model=="TIMEF" || model=="TIM") com.nkappa = 3;
	else if (model=="TVMEF" || model=="TVM") com.nkappa = 4;
	else if (model=="SYM" || model=="GTR") com.nkappa = 5;
		
	//parameters' number: plus another two parameters (Ka/Ks and t)
	com.np = 2 + com.nkappa;

}

GY94::~GY94() {

}


double GY94::rndu(void) {
	w_rndu = w_rndu*69069+1;
	return ldexp((double)w_rndu, -32);
}


/* x[i]=x0[i] + t*p[i] */
double GY94::fun_ls (double t, double x0[], double p[], double x[], int n) {
	int i;  
	for (i=0; i<n; i++) x[i]=x0[i] + t*p[i]; 
	return( lfun2dSdN(x,n) ); 
}



int GY94::PatternWeight () {

   int *fpatti, lst=com.ls,h,ht,j,k=-1;
   int gap=3;
   int nb=1;
   char *zt[2], b1[2]; /* b[][] data at site h */
   double nc = 65;

   //fpatti=(int*)malloc(lst*sizeof(int));
	fpatti= new int[lst];

   for(h=0; h<lst; h++) fpatti[h]=0;

   for(j=0;j<com.ns;j++) {
      zt[j]=(char*)malloc(lst*sizeof(char));
	  //zt[j]= new char[lst];
	}

   com.npatt=0;   

   for(h=0; h<com.ls; h++) {
	   
	   b1[0]=com.z[0][h];
	   b1[1]=com.z[1][h];
	   
	   //cout<<(char)b1[1]<<endl;

	   for(ht=0; ht<com.npatt; ht++) {
		   if(b1[0]==zt[0][ht] && b1[1]==zt[1][ht]) 
			   break;
	   }
	   if(ht==com.npatt) {
		   zt[0][com.npatt]=b1[0];
		   zt[1][com.npatt]=b1[1];

		   ht=com.npatt++;
	   }	
	   
	   fpatti[ht]++;
   } //for (h)
   

   //com.fpatt=(double*)malloc(com.npatt*sizeof(double));
   
   for(h=0;h<com.npatt;h++) {
   	   com.fpatt[h]=fpatti[h];	   
   }
   
   //free(fpatti);
   delete []fpatti;

   for(j=0; j<com.ns; j++) {
	   com.z[j]=(char*)realloc(zt[j],com.npatt*sizeof(char));
   }
   
   return 0;
}

//Preprocess in preparation for estimation
int GY94::preProcess(const char* seq1, const char* seq2) {
	
	int i;
	
//	com.space=NULL;
	com.kappa=2;	com.omega = 0.4;
	
	setmark_61_64();

	com.ls = strlen(seq1);

	for(i=0, snp=0; i<com.ls; i++)
		if (seq1[i]!=seq2[i]) snp++;	
	
	for(i=0; i<com.ns; i++) {      
      com.z[i]=new char[com.ls+1];	  
	}
	
	strcpy(com.z[0], seq1);
	strcpy(com.z[1], seq2);

	com.ls /=3;
    
	EncodeSeqs();	
	PatternWeight();

	return 0;
}

int GY94::setmark_61_64 (void) {

	int i, *code=GeneticCode[com.icode];
	//int c[3],aa0,aa1, by[3]={16,4,1};
	//double nSilent, nStop, nRepl;
	
	Nsensecodon=0;
	for (i=0; i<64; i++) {
		if (code[i]==-1) {
			FROM64[i]=-1; 
		}
		else {
			FROM61[Nsensecodon] = i;
			FROM64[i] = Nsensecodon++;
		}
	}

	com.ncode=Nsensecodon;
	
	return 0;
}

//Convert T,C,A,G to 0,1,2,3
int GY94::transform(char *z, int ls) {
	
	int i;
	char *p;
	
	for (i=0,p=z; i<ls; i++,p++) { 
		*p=(char)convertChar(*p);
	}
	
	return 1;
}

void GY94::EncodeSeqs(void) {

   int j,h,k, b[3];
  
   //T C A G to 0 1 2 3
   for(j=0; j<com.ns; j++) {
      transform(com.z[j],com.ls*3);
   }

   //encode to codon 0-64
    for(j=0; j<com.ns; j++) {
         for(h=0; h<com.ls; h++) {
            b[0]=com.z[j][h*3]; b[1]=com.z[j][h*3+1]; b[2]=com.z[j][h*3+2];
            k=b[0]*16+b[1]*4+b[2];
            com.z[j][h]=(char)FROM64[k];
         }		 
    }
}

int GY94::GetCodonFreqs(double pi[]) {
	int n=com.ncode, i,j,ic,b[3];
	double fb3x4[12], fb4[4];
	int flag[CODON];

	for(i=0, initArray(flag, 64); i<str1.length(); i+=3) {
		j = getID(str1.substr(i,3));
		if (getAminoAcid(j)!='!') flag[j] = 1;
		
		j = getID(str2.substr(i,3));
		if (getAminoAcid(j)!='!') flag[j] = 1;
	}

	//Whether sequences are long enough
	if (sumArray(flag, CODON)==com.ncode) {
		return 0;
	}

	//Codon frequency is estimated from nucleotide frequency(fb3x4) at the three positions
	for(i=0,initArray(fb3x4,12),initArray(fb4,4); i<n; i++) {
		ic=FROM61[i];  
		b[0]=ic/16;
		b[1]=(ic/4)%4; 
		b[2]=ic%4;
		for(j=0; j<3; j++) { 
			fb3x4[j*4+b[j]] += pi[i];  
			fb4[b[j]] += pi[i]/3.; 
		}
	}
	
	//use nul frequencies to get codon frequency
	for (i=0; i<n; i++) {
		ic=FROM61[i];  b[0]=ic/16; b[1]=(ic/4)%4; b[2]=ic%4;
		pi[i]=fb3x4[b[0]]*fb3x4[4+b[1]]*fb3x4[8+b[2]];
	}
	
	scaleArray(1./sumArray(pi,n), pi, n);
	
	return 0;
}

int GY94::gradientB(int n, double x[], double f0, double g[], double space[], int xmark[]) {
/* f0=fun(x) is always provided.
   xmark=0: central; 1: upper; -1: down
*/
   int i,j;
   double *x0=space, *x1=space+n, eh0=Small_Diff, eh;  /* eh0=1e-6 || 1e-7 */

   for(i=0; i<n; i++) {
      eh=eh0*(fabs(x[i])+1);
      if (xmark[i]==0 && SIZEp<1) {    //central 
         for (j=0; j<n; j++)  x0[j]=x1[j]=x[j];
         eh=pow(eh,.67);  x0[i]-=eh; x1[i]+=eh;
         g[i]=(lfun2dSdN(x1,n)-(lfun2dSdN)(x0,n))/(eh*2.0);
      }
      else  {//forward or backward
         for(j=0; j<n; j++)  x1[j]=x[j];
         if (xmark[i]) eh*=-xmark[i];
         x1[i] += eh;
         g[i]=(lfun2dSdN(x1,n)-f0)/eh;
      }
   }
   return(0);
}


double GY94::LineSearch2 (double *f, double x0[], double p[], double step, double limit, double e, double space[], int n) {
/* linear search using quadratic interpolation 
   from x0[] in the direction of p[],
                x = x0 + a*p        a ~(0,limit)
   returns (a).    *f: f(x0) for input and f(x) for output

   x0[n] x[n] p[n] space[n]

   adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
   optimization: An introduction.  Van Nostrand Reinhold Company, New York.
   pp. 62-73.
   step is used to find the bracket and is increased or reduced as necessary, 
   and is not terribly important.
*/
   int ii=0, maxround=10, status, nsymb=0;
   double *x=space, factor=4, small=1e-6, smallgapa=0.2;
   double a0,a1,a2,a3,a4=-1,a5,a6, f0,f1,f2,f3,f4=-1,f5,f6;

   a0=a1=0; f1=f0=*f;
   a2=a0+step; 
   f2=fun_ls(a2, x0,p,x,n);
   if (f2>f1) {
      for (; ;) {
         step/=factor;
         if (step<small) return (0);
         a3=a2;    f3=f2;
         a2=a0+step;  f2=fun_ls(a2, x0,p,x,n);
         if (f2<=f1) break;
      }
   }
   else {       
      for (; ;) {
         step*=factor;
         if (step>limit) step=limit;
         a3=a0+step;  f3=fun_ls(a3, x0,p,x,n);

		 //obtain a bracket
         if (f3>=f2) {	//a1<a2<a3 and f1>f2<f3
			 break;
		 }

         a1=a2; f1=f2;    a2=a3; f2=f3;
         if (step>=limit) {
            *f=f3; return(a3);
         }
      }
   }

   /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
   for (ii=0; ii<maxround; ii++) {	//a4 is the minimum from the parabola over (a1,a2,a3)
	   
	   //2.2.4
	   a4 = (a2-a3)*f1+(a3-a1)*f2+(a1-a2)*f3;	
	   if(fabs(a4)>1e-100) {
		   //2.2.3
		   a4 = ((a2*a2-a3*a3)*f1+(a3*a3-a1*a1)*f2+(a1*a1-a2*a2)*f3)/(2*a4);	
	   }
	   if (a4>a3 || a4<a1) {//out of range: whether a1<a4<a3
		   a4=(a1+a2)/2;		   
	   }
	   else {
		   if((a4<=a2 && a2-a4>smallgapa*(a2-a1)) || (a4>a2 && a4-a2>smallgapa*(a3-a2)))
			   status='Y';
		   else 
			   status='C';
	   }
	   f4 = fun_ls(a4, x0,p,x,n);
	   
	   if (fabs(f2-f4)<e*(1+fabs(f2))) {
		   break;
	   }
	   
	   if (a4<=a2) {    /* fig 2.2.10 */
		   if (a2-a4>smallgapa*(a2-a1)) {
			   if (f4<=f2) { a3=a2; a2=a4;  f3=f2; f2=f4; }
			   else        { a1=a4; f1=f4; }
		   }
		   else {
			   if (f4>f2) {
				   a5=(a2+a3)/2; f5=fun_ls(a5, x0,p,x,n);
				   if (f5>f2) { a1=a4; a3=a5;  f1=f4; f3=f5; }
				   else       { a1=a2; a2=a5;  f1=f2; f2=f5; }
			   }
			   else {
				   a5=(a1+a4)/2; f5=fun_ls(a5, x0,p,x,n);
				   if (f5>=f4)
				   { a3=a2; a2=a4; a1=a5;  f3=f2; f2=f4; f1=f5; }
				   else {
					   a6=(a1+a5)/2; f6=fun_ls(a6, x0,p,x,n);
					   if (f6>f5)
                       { a1=a6; a2=a5; a3=a4;  f1=f6; f2=f5; f3=f4; }
					   else { a2=a6; a3=a5; f2=f6; f3=f5; }
				   }
			   }
		   }
	   }
	   else {                     /* fig 2.2.9 */
		   if (a4-a2>smallgapa*(a3-a2)) {
			   if (f2>=f4) { a1=a2; a2=a4;  f1=f2; f2=f4; }
			   else        { a3=a4; f3=f4; }
		   }
		   else {
			   if (f4>f2) {
				   a5=(a1+a2)/2; f5=fun_ls(a5, x0,p,x,n);
				   if (f5>f2) { a1=a5; a3=a4;  f1=f5; f3=f4; }
				   else       { a3=a2; a2=a5;  f3=f2; f2=f5; }
			   }
			   else {
				   a5=(a3+a4)/2; f5=fun_ls(a5, x0,p,x,n);
				   if (f5>=f4)
				   { a1=a2; a2=a4; a3=a5;  f1=f2; f2=f4; f3=f5; }
				   else {
					   a6=(a3+a5)/2; f6=fun_ls(a6, x0,p,x,n);
					   if (f6>f5)
					   { a1=a4; a2=a5; a3=a6;  f1=f4; f2=f5; f3=f6; }
					   else { a1=a5; a2=a6;  f1=f5; f2=f6; }
				   }
			   }
		   }
	   }
   }
   
   if (f2>f0 && f4>f0)  a4=0;
   if (f2<=f4)  { *f=f2; a4=a2; }
   else         *f=f4;
   
   return a4;
}

double GY94::distance(double x[], double y[], int n) {
	
	int i; 
	double t=0;
	
	for (i=0; i<n; i++) t += square(x[i]-y[i]);
	
	return sqrt(t);
}

int GY94::H_end (double x0[], double x1[], double f0, double f1, double e1, double e2, int n) {
/*   Himmelblau termination rule.   return 1 for stop, 0 otherwise.
*/	
	double r;
	
	if((r=norm(x0,n))<e2)  r=1;
	r*=e1;
	if(distance(x1,x0,n)>=r) 
		return 0;

	r=fabs(f0); 
	if(r<e2) r=1;     
	r*=e1;
	if (fabs(f1-f0)>=r) 
		return 0;

	return 1;
}

int GY94::ming2 (double *f, double x[], double xb[][2], double space[], double e, int n) {
/* n-variate minimization with bounds using the BFGS algorithm
	g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
	xmark[n],ix[n]
	Size of space should be (check carefully?)
	#define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))
	nfree: # free variables
	xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
	x[] has initial values at input and returns the estimates in return.
	ix[i] specifies the i-th free parameter
*/
	int i,j, i1,i2,it, maxround=2000, fail=0, *xmark, *ix, nfree;
	int Ngoodtimes=2, goodtimes=0;
	double small=1e-6, sizep0=0;
	double f0, *g0, *g, *p, *x0, *y, *s, *z, *H, *C, *tv;
	double w,v, alpha, am, h, maxstep=8;
	
	if(n==0) return 0;
	
	g0=space;   g=g0+n;  p=g+n;   x0=p+n;
	y=x0+n;     s=y+n;   z=s+n;   H=z+n;  C=H+n*n; tv=C+n*n;

	xmark = (int*)(tv+2*n);  
	ix = xmark+n;
	
	for(i=0; i<n; i++)  {
		xmark[i]=0;
		ix[i]=i; 
	}
	
	for(i=0,nfree=0;i<n;i++) {
		if(x[i]<=xb[i][0]) {
			x[i]=xb[i][0]; 
			xmark[i]=-1; 
			continue;
		}
		if(x[i]>=xb[i][1]) { 
			x[i]=xb[i][1]; 
			xmark[i]= 1; 
			continue;
		}

		ix[nfree++] = i;
	}
	
	f0=*f=lfun2dSdN(x,n);
	
	copyArray(x,x0,n);
	SIZEp=999; 
	
	gradientB (n, x0, f0, g0, tv, xmark);
	
	initIdentityMatrix(H,nfree);
	
	for(Iround=0; Iround<maxround; Iround++) {
		for(i=0,initArray(p,n); i<nfree; i++)  {
			for(j=0; j<nfree; j++) p[ix[i]] -= H[i*nfree+j]*g0[ix[j]];
		}
		
		sizep0 = SIZEp;
		SIZEp=norm(p,n);      /* check this */
		
		for (i=0,am=maxstep; i<n; i++) {  /* max step length */
			if (p[i]>0 && (xb[i][1]-x0[i])/p[i]<am)
				am=(xb[i][1]-x0[i])/p[i];
			else if (p[i]<0 && (xb[i][0]-x0[i])/p[i]<am)
				am=(xb[i][0]-x0[i])/p[i];
		}
		
		if (Iround==0) {
			h=fabs(2*f0*.01/innerp(g0,p,n));  /* check this?? */
			h=min2(h,am/2000);
		}
		else {
			h=norm(s,nfree)/SIZEp;
			h=max2(h,am/500);
		}
		h=max2(h,1e-5);   h=min2(h,am/5);
		*f=f0;
		alpha = LineSearch2(f,x0,p,h,am, min2(1e-3,e), tv,n); /* n or nfree? */
		
		fail=0;
		for(i=0; i<n; i++)  x[i]=x0[i]+alpha*p[i];
		
		w=min2(2,e*1000);
		if(e<1e-4 && e>1e-6) 
			w=0.01;
		
		if(Iround==0 || SIZEp<sizep0 || (SIZEp<.001 && sizep0<.001)) 
			goodtimes++;
		else
			goodtimes=0;
		if((n==1||goodtimes>=Ngoodtimes) && SIZEp<(e>1e-5?.05:.001) && (f0- *f<0.001) && H_end(x0,x,f0,*f,e,e,n))
			break;		
		
		gradientB (n, x, *f, g, tv, xmark);
		
		/* modify the working set */
		for (i=0; i<n; i++) {         /* add constraints, reduce H */
			
			if (xmark[i]) 
				continue;
			
			if (fabs(x[i]-xb[i][0])<1e-6 && -g[i]<0)  
				xmark[i]=-1;
			else if (fabs(x[i]-xb[i][1])<1e-6 && -g[i]>0)  
				xmark[i]=1;
			
			if (xmark[i]==0) 
				continue;

			copyArray (H, C, nfree*nfree);
			for (it=0; it<nfree; it++) 
				if (ix[it]==i)
					break;
				for (i1=it; i1<nfree-1; i1++) 
					ix[i1]=ix[i1+1];
				for (i1=0,nfree--; i1<nfree; i1++) 
					for (i2=0; i2<nfree; i2++)
					H[i1*nfree+i2]=C[(i1+(i1>=it))*(nfree+1) + i2+(i2>=it)];
		}
		
		for (i=0,it=0,w=0; i<n; i++) {  /* delete a constraint, enlarge H */
			if (xmark[i]==-1 && -g[i]>w)     { 
				it=i; w=-g[i]; 
			}
			else if (xmark[i]==1 && -g[i]<-w) { 
				it=i; w=g[i]; 
			}
		}
		
		if (w>10*SIZEp/nfree) {
			copyArray (H, C, nfree*nfree);

			for (i1=0; i1<nfree; i1++) {
				for (i2=0; i2<nfree; i2++) {
					H[i1*(nfree+1)+i2]=C[i1*nfree+i2];
				}
			}

			for(i1=0; i1<nfree+1; i1++) {
				H[i1*(nfree+1)+nfree]=H[nfree*(nfree+1)+i1]=0;
			}

			H[(nfree+1)*(nfree+1)-1]=1;
			xmark[it]=0;  
			ix[nfree++]=it;
		}
		
		for(i=0,f0=*f; i<nfree; i++) {
			y[i]=g[ix[i]]-g0[ix[i]];
			s[i]=x[ix[i]]-x0[ix[i]]; 
		}
		for (i=0; i<n; i++) {
			g0[i]=g[i];
			x0[i]=x[i]; 
		}
		
		//reduce loop
		for (i=0,w=v=0.; i<nfree; i++) {
			for (j=0,z[i]=0.; j<nfree; j++) {
				z[i]+=H[i*nfree+j]*y[j];
			}
			w+=y[i]*z[i];
			v+=y[i]*s[i];
		}
		if (fabs(v)<small)   { 
			initIdentityMatrix(H,nfree); 
			fail=1; 
			continue; 
		}
		for(i=0; i<nfree; i++)  {
			for(j=0; j<nfree; j++) {
				H[i*nfree+j] += ((1+w/v)*s[i]*s[j]-z[i]*s[j]-s[i]*z[j])/v;
			}
		}
			
	}//end of for(Iround...)
	
	*f=lfun2dSdN(x,n);  
	
	if(nfree==n) { 
		copyArray(H, space, n*n);
		return 1;
	}
   return 0;
}


void GY94::HouseholderRealSym(double a[], int n, double d[], double e[]) {
	
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      } 
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

int GY94::EigenTridagQLImplicit(double d[], double e[], int n, double z[]) {
	
	int m,j,iter,niter=30, status=0, i,k;
	double s,r,p,g,f,dd,c,b, aa,bb;
	
	for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
	for (j=0;j<n;j++) {
		iter=0;
		do {
			for (m=j;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;  /* ??? */
			}
			if (m != j) {
				if (iter++ == niter) {
					status=-1;
					break;
				}
				g=(d[j+1]-d[j])/(2*e[j]);
				
				/* r=pythag(g,1); */
				
				if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
				else                r=sqrt(1+g*g);
				
				g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
				s=c=1;
				p=0;
				for (i=m-1;i>=j;i--) {
					f=s*e[i];
					b=c*e[i];
					
					/*  r=pythag(f,g);  */
					aa=fabs(f); bb=fabs(g);
					if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
					else if(bb==0)             r=0;
					else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }
					
					e[i+1]=r;
					if (r == 0) {
						d[i+1] -= p;
						e[m]=0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<n;k++) {
						f=z[k*n+i+1];
						z[k*n+i+1]=s*z[k*n+i]+c*f;
						z[k*n+i]=c*z[k*n+i]-s*f;
					}
				}
				if (r == 0 && i >= j) continue;
				d[j]-=p; e[j]=g; e[m]=0;
			}
		} while (m != j);
	}
	return status;
}

void GY94::EigenSort(double d[], double U[], int n) {
	
	int k,j,i;
	double p;
	
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i+1;j<n;j++) {
			if (d[j] >= p) p=d[k=j];
		}
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=U[j*n+i];
				U[j*n+i]=U[j*n+k];
				U[j*n+k]=p;
			}
		}
	}
}

int GY94::eigenRealSym(double A[], int n, double Root[], double work[]) {
	
   int status=0;
   HouseholderRealSym(A, n, Root, work);
   status=EigenTridagQLImplicit(Root, work, n, A);
   EigenSort(Root, A, n);

   return(status);
}


int GY94::eigenQREV (double Q[], double pi[], int n, double Root[], double U[], double V[], double spacesqrtpi[]) {

   int i,j, inew, jnew, nnew, status;
   double *pi_sqrt=spacesqrtpi, small=1e-6;

   for(j=0,nnew=0; j<n; j++)
      if(pi[j]>small)
         pi_sqrt[nnew++]=sqrt(pi[j]);

   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

   if(nnew==n) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);

      status=eigenRealSym(U, n, Root, V);
      for(i=0;i<n;i++) for(j=0;j<n;j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0;i<n;i++) for(j=0;j<n;j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i]>small) {
            for(j=0,jnew=0; j<i; j++) 
               if(pi[j]>small) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew] 
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);

      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i]>small ? Root[inew--] : 0);
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>small) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else 
                  V[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>small) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else 
                  U[i*n+j] = (i==j);
            inew--;
         }
         else 
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }

   return status;
}

/********************************************
* Function: parseSubRates
* Input Parameter: string, array of double
* Output: Parse substitution rates according to the given model
* Return Value: int
*********************************************/
int GY94::parseSubRates(string model, double kappa[]) {
	
	/*
	JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
	K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
	TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
	K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
	TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
	TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
	SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
	*/

	kappa[5]=1.0;		//Substitution rate between C and A

	if(model=="JC" || model=="F81") {				
		//Q[i*n+j] = kappa[0] = 1;
		kappa[0]=kappa[1]=kappa[2]=kappa[3]=kappa[4]=kappa[5];
	}
	else if(model=="K2P" || model=="HKY") {//one para.
		//if ((b1+b2)==1 || (b1+b2)==5) Q[i*n+j] = kappa[0];
		kappa[1]=kappa[0];
		kappa[2]=kappa[3]=kappa[4]=kappa[5];
	}
	else if(model=="TNEF" || model=="TN") {//two para.
		//if ((b1+b2)==1)  Q[i*n+j] = kappa[0];
		//else if ((b1+b2)==5)  Q[i*n+j] = kappa[1];
		kappa[2]=kappa[3]=kappa[4]=kappa[5];
	}
	else if(model=="K3P" || model=="K3PUF") {//two para.
		//if ((b1+b2)==1 || (b1+b2)==5)  Q[i*n+j] = kappa[0];
		//else if ((b1+b2)==2 || (b1+b2)==4) Q[i*n+j] = kappa[1];
		kappa[4]=kappa[5];
		kappa[2]=kappa[3]=kappa[1];
		kappa[1] = kappa[0];
	}
	else if(model=="TIMEF" || model=="TIM") {//three para.
		//if ((b1+b2)==1) Q[i*n+j] = kappa[0];
		//else if ((b1+b2)==5)  Q[i*n+j] = kappa[1];
		//else if ((b1+b2)==2 || (b1+b2)==4) Q[i*n+j] = kappa[2];
		kappa[4]=kappa[5];
		kappa[3]=kappa[2];
	}
	else if(model=="TVMEF" || model=="TVM") {//four para.
		//if ((b1+b2)==1 || (b1+b2)==5)  Q[i*n+j] = kappa[0];
		//else if ((b1+b2)==2) Q[i*n+j] = kappa[1];
		//else if ((b1+b2)==4) Q[i*n+j] = kappa[2];
		//else if ((b1==0) && (b2==3)) Q[i*n+j] = kappa[3];
		
		//TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
		//SYM, GTR:   rTC!=rAG!=rTA!=rCG!=rTG!=rCA
		kappa[4]=kappa[3];
		kappa[3]=kappa[2];
		kappa[2]=kappa[1];
		kappa[1]=kappa[0];
	}
	else if(model=="SYM" || model=="GTR") {//five para.
	}
	else {
		return 0;
	}

	return 1;
}

/* Construct Q: transition probability matrix 64*64 */
int GY94::EigenQc (int getstats, double blength, double *S, double *dS, double *dN, double Root[], double U[], double V[], double kappa[], double omega, double Q[]) {
/* This contructs the rate matrix Q for codon substitution and get the eigen
   values and vectors if getstats==0, or get statistics (dS & dN) if 
   getstats==1.  
   The routine is also called by Qcodon2aa for mechanistic amino acid 
   substitution models.
   Input parameters are kappa, omega and com.pi (or com.fb61).

   Statistics calculated include S, dS & dN.
   c0[0,1,2] and c[0,1,2] are rates for the 3 codon positions before and after 
   selection.  c4 is for 4-fold rates.  ts[3] and tv[3] are transition/
   transversion rates, not calculated.

   Under NSsites or other site-class models, this function does not scale 
   Q but calculates the Qfactor_NS.
   DetailOutput() uses this function to calculate dS & dN; for each omega 
   component, dS and dN are not scaled here but are scaled in DetailOutput()
   after Qfactor_NS is calculated.

   aaDist=FIT1 & FIT2:  ap,p*,av,v*, (and w0 for FIT2)
   The argument omega is used only if the model assumes one omega.  For 
   AAClasses, com.pomega is used instead.
*/
	int n=Nsensecodon, i,j,k, ic1,ic2,aa1,aa2;
	int ndiff,pos=0,from[3],to[3];
	double mr, rs0,ra0,rs,ra; /* rho's */
	double d4=0, d0[3],d[3],ts[3],tv[3];  /* rates at positions and 4-fold sites */
	double *pi=com.pi, w=-1, pijQij;
	double space_pisqrt[CODON];
	
	for(i=0; i<3; d[i]=d0[i]=ts[i]=tv[i]=0, i++);
	initArray(Q, n*n);

	//Equal codon frequency
	if(model=="JC" || model=="K2P" || model=="TNEF" || model=="K3P"||
		model=="TIMEF" ||model=="TVMEF"|| model=="SYM") {
		initArray(com.pi, 64, 1.0/(64-getNumNonsense(genetic_code)));
	}

	//Parse substitution rates according to the given model
	parseSubRates(model, kappa);
	
	//Construct Q: transition probability matrix 64*64
	for (i=0, rs0=ra0=rs=ra=0; i<n; i++) {
		
		//codon i
		ic1=FROM61[i]; from[0]=ic1/16; from[1]=(ic1/4)%4; from[2]=ic1%4;
		
		//codon j
		for(j=0; j<i; j++) {
			ic2=FROM61[j]; to[0]=ic2/16; to[1]=(ic2/4)%4; to[2]=ic2%4;
			for(k=0,ndiff=0; k<3; k++) {
				if(from[k]!=to[k]) { 
					ndiff++; pos=k; 
				}
			}

			//consider only one difference between two codons
			if(ndiff!=1) continue;

			int b1 = min2(from[pos],to[pos]);
			int b2 = max2(from[pos],to[pos]);

			/*           01   23   02   13   03   12
			JC, F81:    rTC==rAG =rTA==rCG==rTG==rCA
			K2P, HKY:   rTC==rAG!=rTA==rCG==rTG==rCA
			TNEF,TN:    rTC!=rAG!=rTA==rCG==rTG==rCA
			K3P, K3PUF: rTC==rAG!=rTA==rCG!=rTG==rCA
			TIMEF, TIM:	rTC!=rAG!=rTA==rCG!=rTG==rCA
			TVMEF, TVM: rTC==rAG!=rTA!=rCG!=rTG!=rCA
			SYM, GTR: 	rTC!=rAG!=rTA!=rCG!=rTG!=rCA
			*/
			if     (b1==0 && b2==1) Q[i*n+j]=kappa[0]; /* TC */
			else if(b1==2 && b2==3) Q[i*n+j]=kappa[1]; /* AG */
			else if(b1==0 && b2==2) Q[i*n+j]=kappa[2]; /* TA */
			else if(b1==1 && b2==3) Q[i*n+j]=kappa[3]; /* CG */
			else if(b1==0 && b2==3) Q[i*n+j]=kappa[4]; /* TG */
			else if(b1==1 && b2==2) Q[i*n+j]=kappa[5]; /* CA=1 */

			Q[j*n+i]=Q[i*n+j];
			Q[i*n+j]*=com.pi[j];
			Q[j*n+i]*=com.pi[i]; 
			
			//probability
			pijQij=2*pi[i]*Q[i*n+j];

			aa1=GeneticCode[com.icode][ic1];  
			aa2=GeneticCode[com.icode][ic2];
			if(aa1==aa2) {//synonymous
				rs+=pijQij;
			}
			else {//nonsynonymous
				ra0+=pijQij;
				w = omega;
				
				Q[i*n+j]*=w; Q[j*n+i]*=w;
				ra+=pijQij*w;
			}

		} /* for (j) */
	}    /* for (i) */
	
	mr=rs+ra;
	if(getstats) {
		rs0=rs;
		w=(rs0+ra0);  rs0/=w;  ra0/=w;   *S=rs0*3*com.ls;
		if(blength>=0) {  /* calculates dS & dN */
			if(blength==0) *dS = *dN = 0;
			rs/=mr;
			ra/=mr;
			*dS=blength*rs/(3*rs0);
			*dN=blength*ra/(3*ra0);
		}
		
	}
	else {
		for (i=0; i<n; i++) Q[i*n+i]=-sumArray(Q+i*n,n);
		
		eigenQREV(Q, com.pi, n, Root, U, V, space_pisqrt);
		
		for(i=0; i<n; i++) Root[i]/=mr;
		
	}
	return 0;
}

/* Return maximum-likelihood score */
double GY94::lfun2dSdN(double x[], int np) {
/* likelihood function for calculating dS and dN between 2 sequences,
   com.z[0] & com.z[1]:
         prob(i,j) = PI_i * p(i,j,t)
   
   Data are clean and coded.
   Transition probability pijt is calculated for observed patterns only.
*/
	int n=com.ncode, h,k, ik, z0,z1;
	double  fh,expt[CODON], lnL1=0;
	double *pkappa=com.KAPPA;

	//cout<<name.c_str()<<": "<<com.ncode<<"\t"<<Nsensecodon<<endl;
	
	k=1, ik=0;
	com.kappa=x[k]; 
	for (ik=0; ik<com.nkappa; ik++) 
		pkappa[ik]=x[k++];

	com.omega=x[1+com.nkappa];

	EigenQc(0,-1,NULL,NULL,NULL, Root, U, V, pkappa, com.omega, PMat);
	
	//t = x[0],  exp(Qt)
	for(k=0; k<n; k++) {
		expt[k] = exp(x[0]*Root[k]);
	}

	//com.npatt = number of patterns
	for (h=0; h<com.npatt; h++) {

		if(com.fpatt[h]<Small_Diff) {
			continue;
		}
		
		z0=com.z[0][h]; 
		z1=com.z[1][h];
		
		for(k=0,fh=0;k<n;k++) {
			fh += U[z0*n+k]*expt[k]*V[k*n+z1];
		}
		
		fh*=com.pi[z0];

		lnL1-=log(fh)*com.fpatt[h];
	}

	return lnL1;
}

/* Main function to calculate Ka and Ks, called by "Run" */
int GY94::PairwiseCodon(double space[]) {

	char *pz0[2];
	int npatt0=com.npatt;
	double *fpatt0, ls0=com.ls;
	double fp[CODON*CODON];
	int n=com.ncode, j,k,h;
	double x[10]={.3,1,.5,.5,.5,.5,.3}, xb[10][2]={{1e-6,3.}};
	double kappab[2]={.01,30}, omegab[2]={0.001, 50.};
	double e=1e-6, dS,dN;
	double *pkappa=com.KAPPA;
	
	//fpatt0=(double*)malloc(npatt0*3*sizeof(double));
	fpatt0= new double[npatt0*3];
	for(k=0; k<com.ns; k++) pz0[k]=com.z[k];
	com.z[0]=(char*)(fpatt0+npatt0);  
	com.z[1]=com.z[0]+npatt0;

	copyArray(com.fpatt, fpatt0, npatt0);
	
	//t > snp/length
	if ((snp/length) > 1e-6) xb[0][0] = 3*snp/length;

	//com.nkappa = 2;
	for(j=0; j<com.nkappa; j++) { 
		xb[1+j][0]=kappab[0]; 
		xb[1+j][1]=kappab[1];
	}
	
	k = 1+com.nkappa; 
	xb[k][0] = omegab[0]; 
	xb[k][1] = omegab[1];

	//fp: codon -> codon   frequency
	initArray(fp, CODON*CODON);		
	for(h=0; h<npatt0; h++) {
		j = max2(pz0[0][h], pz0[1][h]);
		k = min2(pz0[0][h], pz0[1][h]);
		fp[j*n+k] += fpatt0[h];
	}	
	
	//fp -> com.fpatt
	for(j=0,com.npatt=0;j<n;j++) {
		for(k=0;k<j+1;k++)
			if(fp[j*n+k]) {
				com.z[0][com.npatt]=(char)j;
				com.z[1][com.npatt]=(char)k;
				com.fpatt[com.npatt++]=fp[j*n+k];
			}
	}
	
	//com.pi: sense codon's frequencies
	for(j=0,initArray(com.pi,n); j<com.npatt; j++) {
		com.pi[(int)com.z[0][j]]+=com.fpatt[j]/(2.*com.ls);
		com.pi[(int)com.z[1][j]]+=com.fpatt[j]/(2.*com.ls);				
	}
	
	GetCodonFreqs(com.pi);	
	
	/* initial values and bounds */
	//divergence time t
	x[0]=-1;
	if(x[0]>3) 
		x[0]=1.5+rndu();
	if(x[0]<1e-6)
		x[0]=.5*rndu();

	//kappas
	for(j=0; j<com.nkappa; j++) 
			x[1+j]=.2+.4*rndu();

	//Ka/Ks
	k=1+com.nkappa;
	x[k]=(3*x[k]+0.6*rndu())/4;
	x[k]=max2(x[k],0.01);
	x[k]=min2(x[k],2);
	
	if ((snp/length) > 1e-6) x[0] = 3.0*snp/length;

	ming2(&lnL, x, xb, space, e, com.np);
	
	EigenQc(1, x[0], &S, &dS, &dN, NULL, NULL, NULL, pkappa, com.omega, PMat);

	Ka = dN;
	Ks = dS;
	N = com.ls*3-S;
	Sd = Ks*S;
	Nd = Ka*N;	
	double scale=(Sd+Nd)/snp;
	Sd /= scale;
	Nd /= scale;

	lnL = -lnL;
	//AICc = -2log(lnL) + 2K + 2K(K+1)/(n-K-1), K=parameters' number, n=sample size
	AICc = -2*lnL + 2.*(com.nkappa+2)*(length/3)/((length/3)-(com.nkappa+2)-1.);

	t = x[0]/3;
	//kappa = com.kappa;
	
	parseSubRates(model, pkappa);
	copyArray(pkappa, KAPPA, NUMBER_OF_RATES);

	k=com.np-1;
	
	com.ls=(int)ls0;
	for(k=0; k<com.ns; k++) com.z[k]=pz0[k];  
	com.npatt=npatt0;
	
	for(h=0; h<npatt0; h++) com.fpatt[h]=fpatt0[h];
	//free(fpatt0);
	delete []fpatt0;
	
	return 0;
}

void GY94::FreeMemPUVR(void){   
   //free(PMat); 
}

/* Main fuction for GY method */
string GY94::Run(const char *seq1, const char *seq2) {
	
	int i;
	
	str1 = seq1;
	str2 = seq2;

	preProcess(seq1, seq2);
	
	com.sspace = max2(800000, 3*com.ncode*com.ncode*(int)sizeof(double));	
	//com.space = (double*)realloc(com.space, com.sspace);
	com.space = new double[com.sspace];

	PairwiseCodon(com.space); 
	
	FreeMemPUVR();
	
	for(i=0; i<com.ns; i++) delete []com.z[i];
	
	//free(com.fpatt);
	//free(com.space);
	delete []com.space;
		
	return parseOutput();
}

