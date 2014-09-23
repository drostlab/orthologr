/*********************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
*
* Filename: MYN.cpp
* Abstract: Definition of Modified YN00 (MYN) class.
*
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Dec.30, 2005

  References: 
	Zhang Zhang, Jun Li, Jun Yu. (2006) Computing Ka and Ks 
	with a consideration of unequal transitional substitutions. 
	BMC Evolutionary Biology, 6:44.
**********************************************************/

#include "MYN.h"

MYN::MYN() {
	name = "MYN";
}

/* Get the two kappas between purines and between pyrimidines */
int MYN::GetKappa(const string seq1, const string seq2) {

	int i,j,k,h,pos,c[2],aa[2],b[2][3],nondeg,fourdeg,by[3]={16,4,1};
	double kappatc_TN[2], kappaag_TN[2], kappa_TN[2];
	double F[2][XSIZE], S[2], wk[2], pi4[4]; 
	double T1, T2, V;//proportions of transitional differences between purines and between
	double kdefault=2, nullValue=NULL;
		
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
		for(j=0; j<16; j++) {
			F[k][j]/=S[k];			
		}
		
		//Transitions between purines
		T1 = 2*F[k][2*DNASIZE+3];
		//Transitions between pyrimidines
		T2 = 2*F[k][0*DNASIZE+1];
		//Tranversions
		V = 1 - T1 - T2 - (F[k][0*DNASIZE+0]+F[k][1*DNASIZE+1]+F[k][2*DNASIZE+2]+F[k][3*DNASIZE+3]);		
		
		//pi[]: the sum probabilty of T, C, A, G, respectively
		for(j=0; j<4; j++) {
			pi4[j] = sumArray(F[k]+j*4,4);
		}
		
		CorrectKappaTN93(S[k], T1, T2, V, pi4, kappatc_TN[k], kappaag_TN[k]);
		wk[k]=((kappatc_TN[k]>0 && kappaag_TN[k]>0)?S[k]:0);
		
		//R = (¦ÐT¦ÐC¦Ê1 + ¦ÐA¦ÐG¦Ê2)/(¦ÐY¦ÐR), kappa = 2R in PAML's DOC
		kappa_TN[k] = 2*(kappatc_TN[k]*pi4[0]*pi4[1] + kappaag_TN[k]*pi4[2]*pi4[3])/((pi4[0]+pi4[1])*(pi4[2]+pi4[3]));

	}

	if(wk[0]+wk[1]==0) {
		kappatc = kappaag = kappa = kdefault;	
	}
	else {
		kappatc = (kappatc_TN[0]*wk[0] + kappatc_TN[1]*wk[1])/(wk[0] + wk[1]);
		kappaag = (kappaag_TN[0]*wk[0] + kappaag_TN[1]*wk[1])/(wk[0] + wk[1]);		
		kappa   = (kappa_TN[0]*wk[0] + kappa_TN[1]*wk[1])/(wk[0] + wk[1]);
	}

	KAPPA[0] = kappatc;
	KAPPA[1] = kappaag;
	
	return 0;
}


/* Calculate transition probability matrix(64*64) */
int MYN::GetPMatCodon(double PMatrix[], double kappa, double omega) {

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
				if ((from[pos]+to[pos]-1)*(from[pos]+to[pos]-5)==0) {
					if (from[pos]+to[pos]==1)  
						PMatrix[i*CODON+j]*=kappatc;//T<->C
					else 
						PMatrix[i*CODON+j]*=kappaag;//A<->G
				}
				
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
		PMatrix[i*CODON+i] =- sumArray(PMatrix+i*CODON,CODON);
		//The sum of transition probability of main diagnoal elements
		mr -= pi[i]*PMatrix[i*CODON+i];
	}
	
	//calculate exp(PMatrix*t)
	eigenQREV(PMatrix, pi, pi_sqrt, CODON, npi0, Root, U, V);
	for(i=0; i<CODON; i++) 
		Root[i]/=mr;
	PMatUVRoot(PMatrix,t,64,U,V,Root);	

	return 0;
}

/* Correct kappas */
int MYN::CorrectKappaTN93(double n, double P1, double P2, double Q, double pi4[], double &kappatc_TN93, double &kappaag_TN93) {

	int failTN93;
	double tc,ag, Y,R, a1=0, a2=0,b=0, A,B,C;
	double Qsmall=min(1e-10,0.1/n), default_kappa = 2, maxkappa = 99;

	kappatc_TN93 = kappaag_TN93 = -1;

	failTN93 = 0;

	Y=pi4[0]+pi4[1];
	R=pi4[2]+pi4[3];

	tc=pi4[0]*pi4[1];
	ag=pi4[2]*pi4[3];


	if ((P1+P2+Q)>1 || (P1+P2)<-1e-10 || Q<-1e-10 || fabs(Y+R-1)>1e-8) {
		return 0;
	}

	if(Q<Qsmall) 
		failTN93 =1;
	else if (Y<=0 || R<=0 || (tc<=0 && ag<=0)) 
		failTN93 = 1;
	else {	//TN93 for multiple substitutions
		A=tc/Y+ag/R;	B=tc+ag;	C=Y*R;
		a1 = 1 - R*P1/(2*ag) - Q/(2*R);
		a2 = 1 - Y*P2/(2*tc) - Q/(2*Y);
		b  = 1 - Q/(2*C);
		if(a1<0 || a2<0 || b<0) {
			failTN93 = 1;			
		}
		else {
			a1 = log(a1);
			a2 = log(a2);
			b  = log(b);
			//Kappa
			kappaag_TN93 = (Y*b - a1)/(-R*b);
			kappatc_TN93 = (R*b - a2)/(-Y*b);
		}
	}
	
	if(failTN93) {	//Fail to correct kappa
		kappatc_TN93 = kappaag_TN93 = default_kappa;
	}
	
	return 1;
}

/* Correct Ka and Ks */
int MYN::CorrectKaksTN93(double n, double P1, double P2, double Q, double pi4[], double &kaks, double &SEkaks) {

	int failTN93;
	double tc,ag, Y,R, a1, a2, b, A, B, C;
	double Qsmall=1e-10;

	a1 = a2 = b = failTN93 = 0;

	Y=pi4[0]+pi4[1];
	R=pi4[2]+pi4[3];

	tc=pi4[0]*pi4[1];
	ag=pi4[2]*pi4[3];

	if (P1+P2+Q>1 || fabs(Y+R-1)>Qsmall || Y<=0 || R<=0 || (tc<=0 && ag<=0)) {		
		failTN93 = 1;		
	}	
	else {	//TN93 for multiple substitutions
		A=tc/Y+ag/R; B=tc+ag; C=Y*R;
		a1 = 1 - R*P1/(2*ag) - Q/(2*R);
		a2 = 1 - Y*P2/(2*tc) - Q/(2*Y);
		b = 1-Q/(2*C);
		if(a1<0 || a2<0 || b<0) {
			failTN93 = 1;			
		}
		else {
			a1 = log(a1);
			a2 = log(a2);
			b  = log(b);			
			//Ka or Ks
			kaks = (-2*ag*a1/R) + (-2*tc*a2/Y) + (-2*(C-ag*Y/R-tc*R/Y)*b);

			double cc1 = 2*ag*R / (2*ag*R-R*R*P1-ag*Q);
			double cc2 = 2*tc*Y / (2*tc*Y-Y*Y*P2-tc*Q);
			double cc3 = 2*ag*ag / (R*(2*ag*R-R*R*P1-ag*Q));
			cc3 += 2*tc*tc / (Y*(2*tc*Y-Y*Y*P2-tc*Q));
			cc3 += (R*R*(Y*Y-2*tc) + Y*Y*(R*R-2*ag)) / (2*R*R*Y*Y-R*Y*Q);
			SEkaks = (square(cc1)*P1 + square(cc2)*P2 + square(cc3)*Q - square(cc1*P1+cc2*P2+cc3*Q))/n;
		}
	}
	
	if (failTN93==1) {	//Use YN00's correction for Ka, Ks
		DistanceF84(n, P1+P2, Q, pi4, Qsmall, kaks, SEkaks);
	}

	return 1;
}

/* Count differences, considering different transitional pathways between purines and between pyrimidines */
int MYN::CountDiffs(const string seq1, const string seq2, double &Sdts1, double &Sdts2, double &Sdtv,double &Ndts1, double &Ndts2, double &Ndtv,double PMatrix[]) {
	int h,i1,i2,i,k, transi, c[2],ct[2], by[3]={16,4,1};
	char aa[2];
	int dmark[3], step[3], b[2][3], bt1[3], bt2[3];
	int ndiff, npath, nstop, sts1path[6], sts2path[6], stvpath[6],nts1path[6], nts2path[6], ntvpath[6];
	double sts1, sts2, stv, nts1, nts2, ntv; /* syn ts & tv, nonsyn ts & tv for 2 codons */
	double ppath[6], sump,p;

	Sdts1=Sdts2=Sdtv=Ndts1=Ndts2=Ndtv=snp=0;
	for (h=0; h<seq1.length(); h+=3)  {		

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
		sts1=sts2=stv=nts1=nts2=ntv=0;
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
			//transi=(transi==1 || transi==5);
			if (aa[0]==aa[1])  { 
				if (transi==5)
					sts1++;
				else if (transi==1)
					sts2++;
				else
					stv++; 
			}
			else {
				if (transi==5)
					nts1++;
				else if (transi==1)
					nts2++;
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

				sts1path[k]=sts2path[k]=stvpath[k]=nts1path[k]=nts2path[k]=ntvpath[k]=0;  
				
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
										
					//ts & tr when syn & nonsyn in path k
					if(aa[0]==aa[1]) { 
						if(transi==5)
							sts1path[k]++;
						else if(transi==1)
							sts2path[k]++;
						else
							stvpath[k]++; 
					}
					else {
						if(transi==5)
							nts1path[k]++;
						else if(transi==1)
							nts2path[k]++;
						else
							ntvpath[k]++; 
					}

					for(i2=0; i2<3; i2++) 
						bt1[i2]=bt2[i2];
				}
				
			}  /* for(k,npath) */
			if (npath==nstop) {  /* all paths through stop codons */
				if (ndiff==2) { 
					nts1 = 0.25;
					nts2 = 0.25;
					ntv  = 1.5; 
				}
				else {
					nts1 = 0.25;
					nts2 = 0.25;
					ntv  = 2.5; 
				}
			}
			else {
				//sum probabilty of all path
				sump=sumArray(ppath,npath);
				if(sump>1e-20) {					
					for(k=0;k<npath;k++) { //p: the probabilty of path k
						p=ppath[k]/sump;
						sts1+=sts1path[k]*p; sts2+=sts2path[k]*p; stv+=stvpath[k]*p;  
						nts1+=nts1path[k]*p; nts2+=nts2path[k]*p; ntv+=ntvpath[k]*p;
					}
				}
			}			
			
		}//end of if(ndiff)

		Sdts1+=sts1;	Sdts2+=sts2;	Sdtv+=stv;
		Ndts1+=nts1;	Ndts2+=nts2;	Ndtv+=ntv;
	}//end of for(h)
	
   return (0);
}

int MYN::DistanceYN00(const string seq1, const string seq2, double &dS,double &dN, double &SEdS, double &SEdN) {

	int j,ir,nround=100, status=1;
	double fbS[4], fbN[4], fbSt[4], fbNt[4];
	double St, Nt, Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv;
	double w0=0, S0=0, N0=0, dS0=0, dN0=0, accu=5e-8, minomega=1e-5, maxomega=99;
	double PMatrix[CODON*CODON];
	
	//initial values for t and omega(Ka/Ks)
	t=0.09; 
	omega=.5;	
	S = N = 0.0;
	initArray(fbS, 4);
	initArray(fbN, 4);
	
	//Count sites of sequence 1		
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

	//Iterative loop
	for (ir=0; ir<nround; ir++) {   /* iteration */

		//Get transition probability matrix from one codon to another
		GetPMatCodon(PMatrix,kappa,omega);	

		//Count differences
		CountDiffs(seq1, seq2, Sdts1, Sdts2, Sdtv, Ndts1, Ndts2, Ndtv, PMatrix);

		//Synonymous(Sd) and nonsynonymous(Nd) differences
		Sd = Sdts1 + Sdts2 + Sdtv;
		Nd = Ndts1 + Ndts2 + Ndtv;

		//Seldom happen
		if (Sd>S) {
			Sdts1 *= (S/Sd);
			Sdts2 *= (S/Sd);
			Sdtv  *= (S/Sd);
		}
		if (Nd>N) {
			Ndts1 *= (N/Nd);
			Ndts2 *= (N/Nd);
			Ndtv  *= (N/Nd);
		}

		//Ks
		CorrectKaksTN93(S, Sdts1/S, Sdts2/S, Sdtv/S, fbS, dS, SEdS);
		//Ka
		CorrectKaksTN93(N, Ndts1/N, Ndts2/N, Ndtv/N, fbN, dN, SEdN);
		

		status=-1;
			
		if(dS<1e-9) {
			status=-1; 
			omega=maxomega; 
		}
		else {
			omega= max(minomega, dN/dS);
		}

		t = dS * 3 * S/(S + N) + dN * 3 * N/(S + N);

		if (fabs(dS-dS0)<accu && fabs(dN-dN0)<accu && fabs(omega-w0)<accu)
			break;
		
		dS0=dS;
		dN0=dN; 
		w0=omega;

	} //end of for(ir) */

	if(ir==nround) 
		status=-2;
	
	return status;
}


/* Count the synonymous and nonsynonymous sites of two sequences */
int MYN::CountSites(const string seq, double &Stot, double &Ntot,double fbS[],double fbN[]) {
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
				if (k+b[j]==1 || k+b[j]==5)	{//transition
					if (k+b[j]==1) r*=kappatc;
					else r*=kappaag;	//(k+b[j]==5)
				}
				
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
