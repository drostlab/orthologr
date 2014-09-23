/************************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: YN00.cpp
* Abstract: Defination of YN00 class.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.31, 2005

* Modified Version: 
* Modified Author:
* Modified Date:

  Note: Source codes are adapted from yn00.c in PAML.

  Reference:
  Yang Z, Nielsen R  (2000)  Estimating Synonymous and 
  Nonsynonymous Substitution Rates Under Realistic 
  Evolutionary Models. Mol Biol Evol 17:32-43.
*************************************************************/

#include "YN00.h"


YN00::YN00() {

	name = "YN";
	initArray(f12pos, 12);
	initArray(pi, 64);
	iteration = 1;
}

void YN00::getFreqency(const string seq1, const string seq2) {
	
	int i;
	double fstop=0.0;

	//Get A,C,G,T frequency at three positions
	for(i=0; i<seq1.length(); i++) {
		f12pos[(i%3)*4+convertChar(seq1[i])]++;			
		f12pos[(i%3)*4+convertChar(seq2[i])]++;
	}
	for(i=0; i<CODONFREQ; i++) 
		f12pos[i]/=(seq1.length()+seq2.length())/3;
	
	//Get 64 amino acid probability	
	for(i=0; i<CODON; i++) {		
		pi[i] = f12pos[i/16] * f12pos[4+(i%16)/4] * f12pos[8+i%4];
		if (getAminoAcid(i)=='!') {
			fstop+=pi[i]; 
			pi[i]=0; 
		}
	}	
	//Scale the sum of pi[] to 1
	for(i=0; i<CODON; i++)
		pi[i]/=(1.0-fstop);

	if (fabs(1-sumArray(pi,CODON))>1e-6) 
		cout<<"Warning: error in get codon freqency."<<endl;

	for(i=0,npi0=0; i<CODON; i++)
		if(pi[i]) 
			pi_sqrt[npi0++]=sqrt(pi[i]);
	
	npi0=CODON-npi0;
	
}


/* Estimate kappa using the fourfold degenerate sites at third codon positions and nondegenerate sites */
int YN00::GetKappa(const string seq1, const string seq2) {

	int i,j,k,h,pos,c[2],aa[2],b[2][3],nondeg,fourdeg,by[3]={16,4,1};
	double ka[2], F[2][XSIZE],S[2],wk[2], T,V, pi4[4];
	double kdefault=2, nullValue=NULL, t;
		
	for(k=0; k<2; k++)
		initArray(F[k],16);
	
	//Get Pi[] of A,C,G,T
	for(h=0; h<seq1.length(); h+=3) {
		
		//c[]: amino acid(0--63)
		c[0]=getID(seq1.substr(h,3));
		c[1]=getID(seq2.substr(h,3));
		//aa[ ]: amino acid
		aa[0]=getAminoAcid(c[0]);
		aa[1]=getAminoAcid(c[1]);		
		//b[][]: 0--3
		for(j=0; j<3; j++) {
			b[0][j] = convertChar(seq1[h+j]);
			b[1][j] = convertChar(seq2[h+j]);
		}		
		
		//Find non-degenerate sites
		for(pos=0; pos<3; pos++) {        
			for(k=0,nondeg=0; k<2; k++) {
				for(i=0; i<4; i++) {
					if(i!=b[k][pos]) 						
						if (getAminoAcid(c[k]+(i-b[k][pos])*by[pos])==aa[k]) 
							break;					
				}
				if (i==4) 
					nondeg++;
			}
			//F[0][]: 0-fold
			if(nondeg==2) {
				F[0][b[0][pos]*4+b[1][pos]]+=.5;
				F[0][b[1][pos]*4+b[0][pos]]+=.5;
			}					
		}
		
		//Find 4-fold degenerate sites at 3rd position
		for(k=0,fourdeg=0;k<2;k++) {
			for(j=0,i=c[k]-b[k][2]; j<4; j++) 
				if(j!=b[k][2] && getAminoAcid(i+j)!=aa[k])
					break;
			if(aa[0]==aa[1] && j==4)
				fourdeg++;
		}
		//F[1][]: 4-fold
		if (fourdeg==2) {
			F[1][b[0][2]*4+b[1][2]]+=.5; 
			F[1][b[1][2]*4+b[0][2]]+=.5;
		}
		
	}//end of for(h)	
	
	for(k=0; k<2; k++) {  /* two kinds of sites */
		
		S[k] = sumArray(F[k],16);
		if(S[k]<=0) { 
			wk[k] = 0; 
			continue; 
		}		
		for(j=0; j<16; j++)
			F[k][j]/=S[k];
		
		//Transition
		T = (F[k][0*4+1]+F[k][2*4+3])*2;
		//Tranversion
		V = 1 - T - (F[k][0*4+0]+F[k][1*4+1]+F[k][2*4+2]+F[k][3*4+3]);		
		
		//pi[]: the sum probabilty of T, C, A, G, respectively
		for(j=0; j<4; j++) {
			pi4[j] = sumArray(F[k]+j*4,4);
		}

		//Correct kappa
		DistanceF84(S[k],T,V,pi4, ka[k], t, nullValue);
		wk[k]=(ka[k]>0?S[k]:0);	
	}

	if(wk[0]+wk[1]==0) {
		kappa=kdefault;
	}
	else {
		kappa=(ka[0]*wk[0]+ka[1]*wk[1])/(wk[0]+wk[1]);
	}

	KAPPA[0] = KAPPA[1] = kappa;
	
	return 0;
}


/* Correct for multiple substitutions */
int YN00::DistanceF84(double n, double P, double Q, double pi4[], double &k_HKY, double &t, double &SEt) {

	int failF84=0,failK80=0,failJC69=0;
	double tc,ag, Y,R, a=0,b=0, A,B,C, k_F84;
	double Qsmall=min(1e-10,0.1/n), maxkappa=2,maxt=99;
	
	k_HKY=-1;

	Y=pi4[0]+pi4[1];
	R=pi4[2]+pi4[3];

	tc=pi4[0]*pi4[1];
	ag=pi4[2]*pi4[3];

	if (P+Q>1) {
		t=maxt; 
		k_HKY=1; 
		return 3;
	}
	if (P<-1e-10 || Q<-1e-10 || fabs(Y+R-1)>1e-8) {
		return 3;		
	}

	//HKY85
	if(Q<Qsmall) 
		failF84=failK80=1;
	else if (Y<=0 || R<=0 || (tc<=0 && ag<=0)) 
		failF84=1;
	else {
		A=tc/Y+ag/R; B=tc+ag; C=Y*R;
		a=(2*B+2*(tc*R/Y+ag*Y/R)*(1-Q/(2*C)) - P) / (2*A);
		b=1-Q/(2*C);
		if (a<=0 || b<=0)
			failF84=1;
	}

	if (!failF84) {
		a=-.5*log(a); b=-.5*log(b);
		if(b<=0) 
			failF84=1;
		else {
			k_F84=a/b-1;
			t = 4*b*(tc*(1+ k_F84/Y)+ag*(1+ k_F84/R)+C);
			k_HKY = (B + (tc/Y+ag/R)* k_F84)/B; /* k_F84=>k_HKY85 */
			
			//Standard errors
			a = A*C/(A*C-C*P/2-(A-B)*Q/2);
			b = A*(A-B)/(A*C-C*P/2-(A-B)*Q/2) - (A-B-C)/(C-Q/2);
			SEt = sqrt((a*a*P+b*b*Q-square(a*P+b*Q))/n);			
		}
	}

	//K80
	if(failF84 && !failK80) {  /* try K80 */
		a=1-2*P-Q;  b=1-2*Q;
		if (a<=0 || b<=0)
			failK80=1;
		else {
			a=-log(a); b=-log(b);
			if(b<=0)
				failK80=1;
			else {
				k_HKY=(.5*a-.25*b)/(.25*b);
				t = .5*a+.25*b;
			}
			if(SEt) {
				a=1/(1-2*P-Q); b=(a+1/(1-2*Q))/2;
				SEt = sqrt((a*a*P+b*b*Q-square(a*P+b*Q))/n);
			}
		}
	}

	if(failK80) {/* try JC69 */
		if((P+=Q)>=.75) {
			failJC69=1; 
			P=.75*(n-1.)/n; 
		}
		t = -.75*log(1-P*4/3.); 

		if(t>99)
			t=maxt;
		if(SEt) {
			SEt = sqrt(9*P*(1-P)/n) / (3-4*P);
		}
	}

	if(k_HKY>99) 
		k_HKY=maxkappa;
	
	return(failF84 + failK80 + failJC69);
}

int YN00::DistanceYN00(const string seq1, const string seq2, double &dS,double &dN, double &SEKs, double &SEKa) {

	int j,k,ir,nround=10, status=0;
	double fbS[4], fbN[4], fbSt[4], fbNt[4], St, Nt, Sdts, Sdtv, Ndts, Ndtv, k_HKY;
	double w0=0, dS0=0, dN0=0, accu=5e-4, minomega=1e-5, maxomega=99;
	double PMatrix[CODON*CODON];
	
	if(t==0)
		t=.5;  
	if(omega<=0) 
		omega=1;
	for(k=0; k<4; k++) 
		fbS[k]=fbN[k]=0;

	//Count sites of sequence 1
	S=N=0;
	CountSites(seq1, St, Nt, fbSt, fbNt);
	S+=St/2; 
	N+=Nt/2;
	for(j=0; j<4; j++) {
		fbS[j]+=fbSt[j]/2; 
		fbN[j]+=fbNt[j]/2; 
	}
	//Count sites of sequence 2
	CountSites(seq2, St, Nt, fbSt, fbNt);
	S+=St/2; 
	N+=Nt/2;
	for(j=0; j<4; j++) {
		fbS[j]+=fbSt[j]/2; 
		fbN[j]+=fbNt[j]/2; 
	}
	
	if (t<0.001 || t>5) 
		t=0.5; 
	if (omega<0.01 || omega>5)
		omega=.5;

	
	for (ir=0; ir<(iteration?nround:1); ir++) {   /* iteration */
		if(iteration) 
			GetPMatCodon(PMatrix,kappa,omega);
		else
			for(j=0; j<CODON*CODON; j++) {
				PMatrix[j]=1;
			}
		
		CountDiffs(seq1,seq2, Sdts, Sdtv, Ndts, Ndtv, PMatrix);
		
		Sd = Sdts + Sdtv;
		Nd = Ndts + Ndtv;

		//synonymous
		DistanceF84(S, Sdts/S, Sdtv/S, fbS, k_HKY, dS, SEKs);
		//nonsynonymous
		DistanceF84(N, Ndts/N, Ndtv/N, fbN, k_HKY, dN, SEKa);
			
		if(dS<1e-9) {
			status=-1; 
			omega=maxomega; 
		}
		else {
			omega= max(minomega, dN/dS);
		}

		t = dS * 3 * S/(S + N) + dN * 3 * N/(S + N);
		
		if ( fabs(dS-dS0)<accu && fabs(dN-dN0)<accu && fabs(omega-w0)<accu )
			break;
		
		dS0=dS;
		dN0=dN; 
		w0=omega;

	} //end of for(ir) */

	if(ir==nround) 
		status=-2;
	
	return status;
}

//Count differences between two compared codons
int YN00::CountDiffs(const string seq1, const string seq2, double &Sdts,double &Sdtv,double &Ndts, double &Ndtv,double PMatrix[]) {
	int h,i1,i2,i,k, transi, c[2],ct[2], by[3]={16,4,1};
	char aa[2];
	int dmark[3], step[3], b[2][3], bt1[3], bt2[3];
	int ndiff, npath, nstop, stspath[6],stvpath[6],ntspath[6],ntvpath[6];
	double sts,stv,nts,ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
	double ppath[6], sump,p;

	snp = 0;
	for (h=0,Sdts=Sdtv=Ndts=Ndtv=0; h<seq1.length(); h+=3)  {		

		c[0]=getID(seq1.substr(h,3));
		c[1]=getID(seq2.substr(h,3));
		//Difference?
		if (c[0]==c[1])
			continue;

		for(i=0; i<2; i++) {
			b[i][0]=c[i]/16;
			b[i][1]=(c[i]%16)/4; 
			b[i][2]=c[i]%4;
			aa[i]=getAminoAcid(c[i]);
		}

		//ndiff: differences of two codons
		ndiff=0;
		sts=stv=nts=ntv=0;
		//dmark[]: position of different codon 
		for(k=0; k<3; k++) {
			dmark[k] = -1;
			if (b[0][k]!=b[1][k]) 
				dmark[ndiff++]=k;
		}

		snp+=ndiff;

		npath=1;
		if(ndiff>1)
			npath=(ndiff==2)?2:6;

		if (ndiff==1) {
			transi=b[0][dmark[0]]+b[1][dmark[0]];
			transi=(transi==1 || transi==5);
			if (aa[0]==aa[1])  { 
				if (transi)
					sts++; 
				else
					stv++; 
			}
			else {
				if (transi) 
					nts++;
				else
					ntv++;
			}
		}
		else {   /* ndiff=2 or 3 */
			
			nstop=0;
			for(k=0; k<npath; k++) {
				
				//set the step[]
				for(i1=0; i1<3; i1++)
					step[i1]=-1;
				if (ndiff==2) {
					step[0]=dmark[k];
					step[1]=dmark[1-k];  
				}
				else {
					step[0]=k/2; 
					step[1]=k%2;
					if (step[0]<=step[1])
						step[1]++;
					step[2]=3-step[0]-step[1];
				}//end of set the step[]

				for(i1=0; i1<3; i1++)
					bt1[i1]=bt2[i1]=b[0][i1];

				stspath[k]=stvpath[k]=ntspath[k]=ntvpath[k]=0;  
				
				//ppath[]: probabilty of each path
				for (i1=0,ppath[k]=1; i1<ndiff; i1++) { 
					bt2[step[i1]] = b[1][step[i1]];

					//ct[]: mutated codon's ID(0--63)
					for (i2=0,ct[0]=ct[1]=0; i2<3; i2++) {
						ct[0]+=bt1[i2]*by[i2];
						ct[1]+=bt2[i2]*by[i2];
					}
					//ppath[k]: probabilty of path k
					ppath[k]*=PMatrix[ct[0]*CODON+ct[1]];
					for(i2=0; i2<2; i2++) 
						aa[i2]=getAminoAcid(ct[i2]);
					
					if (aa[1]=='!') {
						nstop++;  
						ppath[k]=0;
						break;
					}

					transi=b[0][step[i1]]+b[1][step[i1]];
					transi=(transi==1 || transi==5);  /* transition? */
					
					//ts & tr when syn & nonsyn in path k
					if(aa[0]==aa[1]) { 
						if(transi)
							stspath[k]++;
						else
							stvpath[k]++; 
					}
					else {
						if(transi)
							ntspath[k]++; 
						else
							ntvpath[k]++; 
					}

					for(i2=0; i2<3; i2++) 
						bt1[i2]=bt2[i2];
				}
				
			}  /* for(k,npath) */
			if (npath==nstop) {  /* all paths through stop codons */
				if (ndiff==2) { 
					nts=.5; 
					ntv=1.5; 
				}
				else {
					nts=.5; 
					ntv=2.5; 
				}
			}
			else {
				//sum probabilty of all path
				sump=sumArray(ppath,npath);
				if(sump>1e-20) {					
					for(k=0;k<npath;k++) { //p: the probabilty of path k
						p=ppath[k]/sump;
						sts+=stspath[k]*p; stv+=stvpath[k]*p;  
						nts+=ntspath[k]*p; ntv+=ntvpath[k]*p;
					}
				}
			}			
			
		}//end of if(ndiff)
		Sdts+=sts;
		Sdtv+=stv;
		Ndts+=nts;
		Ndtv+=ntv;

	}//end of for(h)
	
   return (0);
}

/* Calculate transition probability matrix(64*64) */
int YN00::GetPMatCodon(double PMatrix[], double kappa, double omega)
{
/* Get PMat=exp(Q*t) for weighting pathways
	*/
	int i,j,k, ndiff,pos=0,from[3],to[3];
	double mr;
	char c[2];
	double U[CODON*CODON], V[CODON*CODON], Root[CODON*CODON];
	
	initArray(PMatrix, CODON*CODON);
	initArray(U, CODON*CODON);
	initArray(V, CODON*CODON);
	initArray(Root, CODON*CODON);
	
	for(i=0; i<CODON; i++) {
		for(j=0; j<i; j++) {

			//codon 'from'
			from[0]=i/16; from[1]=(i/4)%4; from[2]=i%4;
			//codon 'to'
			to[0]=j/16;   to[1]=(j/4)%4;   to[2]=j%4;
			//amino acid of 'from' and 'to'
			c[0]=getAminoAcid(i);
			c[1]=getAminoAcid(j);
			//stop codon
			if (c[0]=='!' || c[1]=='!')  
				continue;

			//whether two codons only have one difference
			for (k=0,ndiff=0; k<3; k++) {
				if (from[k]!=to[k]) {
					ndiff++; 
					pos=k; 
				}
			}			
			if (ndiff==1) {
				//only have one difference
				PMatrix[i*CODON+j]=1;
				//transition
				if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0)  
					PMatrix[i*CODON+j]*=kappa;
				
				//nonsynonymous
				if(c[0]!=c[1])
					PMatrix[i*CODON+j]*=omega;
				
				//diagonal element is equal
				PMatrix[j*CODON+i]=PMatrix[i*CODON+j];
			}
		}			
	}

	//PMatrix[](*Q): transition probability matrix
	for(i=0; i<CODON; i++) 
		for(j=0; j<CODON; j++) 
			PMatrix[i*CODON+j]*=pi[j];

	//scale the sum of PMat[][j](j=0-63) to zero
	for (i=0,mr=0; i<CODON; i++) { 
		PMatrix[i*CODON+i]=-sumArray(PMatrix+i*CODON,CODON);
		//The sum of transition probability of main diagnoal elements
		mr-=pi[i]*PMatrix[i*CODON+i];
	}

	//calculate exp(PMatrix*t)
	eigenQREV(PMatrix, pi, pi_sqrt, CODON, npi0, Root, U, V);
	for(i=0; i<CODON; i++) 
		Root[i]/=mr;
	PMatUVRoot(PMatrix, t, CODON, U, V, Root);	

	return 0;
}



/* P(t) = U * exp{Root*t} * V */
int YN00::PMatUVRoot (double P[], double t, int n, double U[], double V[], double Root[]) {

	int i,j,k;
	double expt, uexpt, *pP;
	double smallp = 0;
	
	
	for (k=0,initArray(P,n*n); k<n; k++)
		for (i=0,pP=P,expt=exp(t*Root[k]); i<n; i++)
			for (j=0,uexpt=U[i*n+k]*expt; j<n; j++)
				*pP++ += uexpt*V[k*n+j];
			
	for(i=0;i<n*n;i++)
		if(P[i]<smallp) 
			P[i]=0;							
				
	return (0);
}


/* Count the synonymous and nonsynonymous sites of two sequences */
int YN00::CountSites(const string seq, double &Stot, double &Ntot,double fbS[],double fbN[]) {
	int h,i,j,k, c[2],aa[2], b[3], by[3]={16,4,1};
	double r, S,N;
	
	Stot=Ntot=0;  
	initArray(fbS, 4);
	initArray(fbN, 4);

	for (h=0; h<seq.length(); h+=3) {

		//Get codon id and amino acid
		c[0]=getID(seq.substr(h,3));
		aa[0]=getAminoAcid(c[0]); 		
		for(i=0; i<3; i++) {
			b[i]=convertChar(seq[h+i]); 
		}		

		for (j=0,S=N=0; j<3; j++) {
			for(k=0; k<4; k++) {    /* b[j] changes to k */
				if (k==b[j]) 
					continue;
				//c[0] change at position j
				c[1] = c[0]+(k-b[j])*by[j];
				aa[1] = getAminoAcid(c[1]);

				if(aa[1]=='!') 
					continue;
				
				r=pi[c[1]];				
				if (k+b[j]==1 || k+b[j]==5)	//transition
					r*=kappa;
				
				if (aa[0]==aa[1]) { //synonymous
					S+=r;
					fbS[b[j]]+=r; //syn probability of A,C,G,T					
				}
				else { //nonsynonymous
					N+=r;
					fbN[b[j]]+=r; //nonsyn probability of A,C,G,T					
				}
			}
		}
		Stot+=S;
		Ntot+=N;
	}
	
	//Scale Stot+Ntot to seq.length()
	r=seq.length()/(Stot+Ntot);
	Stot*=r; 
	Ntot*=r;

	//get probablity of syn of four nul.
	r=sumArray(fbS,4);
	for(k=0; k<4; k++)
		fbS[k]/=r;

	//get probablity of nonsyn of four nul.
	r=sumArray(fbN,4);
	for(k=0; k<4; k++)
		fbN[k]/=r;

	return 0;
}


string YN00::Run(string seq1, string seq2) {	

	t=0.4; 
	kappa = NA;
	omega = 1;	
	Ks = Ka = 0.1;
	
	getFreqency(seq1, seq2);
	GetKappa(seq1, seq2);
	DistanceYN00(seq1, seq2, Ks, Ka, SEKs, SEKa);
	
	t = (S*Ks+N*Ka)/(S+N);

	return parseOutput();
}



//The following functions are used to calculate PMatrix by Taylor.

int YN00::eigenQREV (double Q[], double pi[], double pi_sqrt[], int n, int npi0,
               double Root[], double U[], double V[])
{
/* 
This finds the eigen solution of the rate matrix Q for a time-reversible 
Markov process, using the algorithm for a real symmetric matrix.
Rate matrix Q = S * diag{pi} = U * diag{Root} * V, 
where S is symmetrical, all elements of pi are positive, and U*V = I.
pi_sqrt[n-npi0] has to be calculated before calling this routine.

  [U 0] [Q_0 0] [U^-1 0]    [Root  0]
  [0 I] [0   0] [0    I]  = [0     0]
  
	Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
	*/
	int i,j, inew,jnew, nnew=n-npi0, status;
	
	//npi0 is the number of stop codons in selected genetic table

	if(npi0==0) {	//seldom occur
		
		//Set U[64*64]
		for(i=0; i<n; i++) {			
			for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++) 
				U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);
		}
		
		//Set U[64*64]
		status=eigenRealSym(U, n, Root, V);

		for(i=0;i<n;i++) 
			for(j=0;j<n;j++)
				V[i*n+j] = U[j*n+i] * pi_sqrt[j];

		for(i=0;i<n;i++) 
			for(j=0;j<n;j++) 
				U[i*n+j] /= pi_sqrt[i];
	}
	else {
		for(i=0,inew=0; i<n; i++) {
			if(pi[i]) {
				for(j=0,jnew=0; j<i; j++) 
					if(pi[j]) {
						U[inew*nnew+jnew] = U[jnew*nnew+inew] = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
						jnew++;
					}
				U[inew*nnew+inew] = Q[i*n+i];
				inew++;
			}
		}
		status=eigenRealSym(U, nnew, Root, V);
		
		for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
			Root[i] = (pi[i] ? Root[inew--] : 0);

		for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
			if(pi[i]) {
				for(j=n-1,jnew=nnew-1; j>=0; j--)
					if(pi[j]) {
						V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
						jnew--;
					}
					else 
						V[i*n+j] = (i==j);
				inew--;
			}
			else 
				for(j=0; j<n; j++)  
					V[i*n+j] = (i==j);
		}

		for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
			if(pi[i]) {
				for(j=n-1,jnew=nnew-1;j>=0;j--)
					if(pi[j]) {
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
	
	return(status);
}

void YN00::HouseholderRealSym(double a[], int n, double d[], double e[]) {
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix 
a[n*n] into a tridiagonal matrix represented by d and e.
d[] is the diagonal (eigends), and e[] the off-diagonal.
	*/
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



int YN00::EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.  
d[] is the diagonal (eigenvalues), e[] is the off-diagonal
z[n*n]: as input should have the identity matrix to get the eigen solution of the 
tridiagonal matrix, or the output from HouseholderRealSym() to get the 
eigen solution to the original real symmetric matrix.
z[n*n]: has the orthogonal matrix as output

  Adapted from routine tqli in Numerical Recipes in C, with reference to
  LAPACK fortran code.
  Ziheng Yang, May 2001
	*/
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
	return(status);
}

void YN00::EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
	*/
	int k,j,i;
	double p;
	
	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i+1;j<n;j++)
			if (d[j] >= p) p=d[k=j];
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

int YN00::eigenRealSym(double A[], int n, double Root[], double work[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return, 
A has the right vectors and Root has the eigenvalues. work[n] is the working space.
The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(), 
and then using the QL algorithm with implicit shifts.  

  Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
  Ziheng Yang, 23 May 2001
	*/
	int status=0;
	HouseholderRealSym(A, n, Root, work);
	status=EigenTridagQLImplicit(Root, work, n, A);
	EigenSort(Root, A, n);
	
	return(status);
}

