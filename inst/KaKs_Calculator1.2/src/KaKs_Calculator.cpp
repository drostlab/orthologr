/***************************************************************
* Copyright (c) 2006, BGI of Chinese Academy of Sciences
* All rights reserved.
 
* Filename: KaKs_Calculator.cpp
* Abstract: including maximum-likelihood and approximate methods.

* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Apr.2, 2006
****************************************************************/

#include "KaKs.h"

int main(int argc, const char* argv[]) {

	try {
		KAKS kk;
		if(!kk.Run(argc, argv)) throw 1;
	}
	catch (...) {
		
	}
	return 0;
}



