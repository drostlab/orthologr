/************************************************************
* Copyright (C) 2005, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: ConcatenatePairs.cpp
* Abstract: Concatenate all pairs of sequences with AXT format
			to a pair of sequences (maybe very long).

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Feb.2, 2005

* Modified Version: 
* Modified Author:
* Modified Date:
*************************************************************/
#include<string>
#include<iostream>
#include<fstream>

using namespace std;

string result;	//Result for outputing into a file

bool readFile(const char *filename) {

	bool flag = true;	
	string seq1, seq2;
	seq1 = seq2 ="";

	try	{
		ifstream is(filename);
		if (!is) {
			cout<<"\nError in opening file..."<<endl;
			throw 1;
		}
		
		cout<<"\nPlease wait while reading sequences and concatenating..."<<endl;

		string temp="", name="", str="";
		
		while (getline(is, temp, '\n'))	{
			
			name = temp;
			
			getline(is, temp, '\n');
			while (temp!="") {				
				str += temp;
				getline(is, temp, '\n');				
			}

			seq1 += str.substr(0, str.length()/2);
			seq2 += str.substr(str.length()/2, str.length()/2);
			name = str = "";			
		}
		
		is.close();
		is.clear();	
		
		result += "Concatenate-";	result += filename;	result += '\n';	
		result += seq1;		result += '\n';
		result += seq2;		result += '\n';
		result += '\n';

	}
	catch (...) {
		
		flag = false;
	}
	
	return flag;
}

bool writeFile(const char *filename, const char *content) {
	
	bool flag=true;

	try {	
		ofstream os(filename);
		os<<content;
		os.close();		
	}
	catch (...) {
		cout<<"Error in writing file..."<<endl;
		flag = false;
	}
	
	return flag;
}


int main(int argc, char* argv[]) {

	if (argc!=3) {
		cout<<"Description: Concatenate all pairs of sequences with AXT format to a pair of sequences."<<endl;
		cout<<"Usage: ConPairs [AXT Filename] [Output Filename]"<<endl;
		return 0;
	}

	try {
		if (!readFile(argv[1])) throw 1;
		if (!writeFile(argv[2], result.c_str())) throw 1;
		cout<<"Mission accomplished."<<endl;

	}
	catch (...) {
		cout<<"Mission failed."<<endl;
	}

	return 1;
}
