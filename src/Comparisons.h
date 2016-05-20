/*

Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu

Remove the brackets to email me.

This file is part of libsequence.

libsequence is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libsequence is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
long with libsequence.  If not, see <http://www.gnu.org/licenses/>.


Modified by Sarah Scharfenberg and Hajk-Georg Drost 2014 to work 
in orthologr without using external libraries from libsequence.

 All changes are also free under the terms of GNU General Public License
 version 2 of the License, or any later version.

*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <assert.h>
#include <string>

#include "Sequence.h"
   
   
// [[Rcpp::export]]
bool Different (const std::string& seq1, const std::string& seq2, bool skip_missing = false, bool nucleic_acid = false)

{
      if(! (seq1.length () == seq2.length ()))
	return true;
      if (skip_missing)
      {
          char MISSING = 'N';
          if (!nucleic_acid)
            MISSING = 'X';//for peptide sequences
          for (unsigned i = 0 ; i < seq1.length() ; ++i)
          {
	      const char _ch1 = char(toupper(seq1[i])),
		_ch2 = char(toupper(seq2[i]));
	      if( (_ch1 != MISSING && _ch2 != MISSING) && _ch1 != _ch2 )
		return true;
	  }
      }
      else
	{
          for (unsigned i = 0 ; i < seq1.length() ; ++i)
	    {
	      if ( char(toupper(seq1[i])) != char(toupper(seq2[i])) ) 
		return true;
	    }
	}
      return false;
}

// [[Rcpp::export]]
unsigned NumDiffs (const std::string & seq1, const std::string & seq2, bool skip_missing = false, bool nucleic_acid = false)
    
{
      unsigned ndiff = 0;
      size_t len = seq1.length();
      if (seq1.length() != seq2.length())
	{
	  len = (seq1.length() < seq2.length()) ? seq1.length() : seq2.length();
	}
      char MISSING = 'N';
      if (!nucleic_acid)
	MISSING = 'X';//for peptide sequences
      for (int i = 0; unsigned (i) < len; ++i)
        {
	  const char _ch1 = char(toupper(seq1[i])), 
	    _ch2 = char(toupper(seq2[i]));
          if(skip_missing == true)
	    {
	      if( (_ch1 != MISSING && _ch2 != MISSING) && (_ch1 != _ch2 ) )
		++ndiff;
	    }
	  else
	    {
              if (_ch1 != _ch2)
                ++ndiff;
	    }
        }
      return ndiff;
}
    
// [[Rcpp::export]]
char toChar(int x) {
   return x;
}

    
// [[Rcpp::export]]
std::string TsTv (int i, int j)
{
      if( ! (i<=3 && j <= 3)){ 
	// catch a call of TsTv without using nucToInt. 
	i = nucToInt( toChar(i));
	j = nucToInt( toChar(j));
	if( ! (i<=3 && j <= 3)){
		Rcpp::Rcout << "Wrong Nucleotide "<<"i="<<i<<" j="<<j<<endl;
		return "Unknown";
	}
      }
      int type = i + j;
      if (type%2!=0.)        //if odd
        {
          return ("Tv");	//a transversion
        }
      else if (type%2==0.)	//if even
        {
          return ("Ts");	//a transition
        }

      Rcpp::Rcout << "Neither Transition, not Transversion" <<endl;
      return ("Unknown");	//can be used for error checking
}
    
    
// [[Rcpp::export]]
bool NotAGap(char c)
/*
 \return true if a c is not a gap character, false otherwise.
 \note Currently, only '-' is considered to be a gap character
*/
{
      switch(c)
        {
        case '-':
          return false;//false, it is a gap...
          break;
        default:
          return true;
          break;
        }
      return true;
}

