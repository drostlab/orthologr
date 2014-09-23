/************************************************************
* Copyright (c) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: KaKs.h
* Abstract: Declaration of KAKS class including several methods.

* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

*************************************************************/
#if !defined(KAKS_H)
#define  KAKS_H

#include "base.h"
#include "NG86.h"
#include "LWL85.h"
#include "LPB93.h"
#include "GY94.h"
#include "YN00.h"
#include "MYN.h"
#include "MSMA.h"

using namespace std;

/* KAKS class */
class KAKS: public Base {
	
public:	
	KAKS();
	~KAKS();

	/* Main function to call kaks's methods */
	bool Run(int argc, const char* argv[]);		
	/* Read and Calculate seq, called in "Run" main function */
	bool ReadCalculateSeq(string filename);
	
	/* Initialize class, ready for running */
	int Initialize();
	/* Unitialize class, for unloading */
	int Uninitialize();

protected:		
	/* Use several methods to calculate ka/ks */
	bool calculateKaKs();
	/* Show help information */
	void helpInfo();

	/* NONE: an in-house algorithm in BGI, that is NG86 without correction */
	void start_NONE();
	/* NG86 */
	void start_NG86();
	/* LWL85 */
	void start_LWL85();
	/* Modified LWL85 */
	void start_MLWL85();
	/* LPB93 */
	void start_LPB93();
	/* Modified LPB93 */
	void start_MLPB93();
	/* GY94 */
	void start_GY94();	
	/* YN00 */
	void start_YN00();
	/* MYN */
	void start_MYN();	
	/* Model Selection and Model Averaging */
	void start_MSMA();
	

	/* Get GCC of entire sequences and of three codon positions */
	void getGCContent(string str);
	/* Check the sequence whether is valid or not */
	bool checkValid(string name, string str1, string str2);
	/* Parse the input parameters */
	bool parseParameter(int argc, const char* argv[]);
	/* Show input parameters' information on screen */
	void showParaInfo();
	/* Get title information for writing into file */
	string getTitleInfo();

public:
	/* Methods' name and reference */
	vector<string> method_name;
	vector<string> method_ref;
	
	/* Parameters' title in outputing file */
	vector<string> titleInfo;

	/* Results for windows parser that shows results on ListCtrl */
	string result4Win;
	
	/* File name for output */
	string output_filename;
	/* Sequence file name */
	string seq_filename;

	/* Flag for whether to run NG86, MLWL85, MLPB93, GY94, YN00, MYN, MS/A=model selection/averaging */
	bool none, ng86, lwl85, lpb93, yn00, mlwl85, mlpb93, gy94, myn, ms, ma;	
	/* Number of compared pairwise sequences */
	unsigned long number;	//Maybe too many

protected:
	/* File name for detailed results for model selection */
	string detail_filename;
	/* Detailed results */
	string details; 
	
private:
	/* The temporary results for write into file */
	string result;
	/* Output stream */
	ofstream os;
	/* A pair of sequence */
	string seq1, seq2;
}; 

#endif



