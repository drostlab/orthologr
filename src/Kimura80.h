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
version 3 of the License, or any later version.

*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "Comparisons.h"

class Kimura80
{
    private:

      unsigned num_Ts, num_Tv;
      size_t seqlen;        	//total sequence length
      size_t sites_compared;	//number of ungapped sites in the data
      
      double divergence, P, Q;
      
      
      void Compute (const pair<string,string> *seq1, 
                    const pair<string,string> *seq2)
      {
	  unsigned i;
		
	  unsigned ungapped_sites = 0;
	  for (i = 0; i < seqlen; ++i)        //iterate over the sequence
	    {
	      string type = "";
	      if (NotAGap(seq1->second[i]) && NotAGap(seq2->second[i]))
		{
	++ungapped_sites;
			    //if the sites differ at that position
	if (toupper(seq1->second[i]) != toupper(seq2->second[i]))	
	  {
				//check if difference is Ts or Tv
	    type = TsTv (nucToInt(seq1->second[i]), nucToInt(seq2->second[i]));	
	  
				if (type == "Ts")
				  {	//Ts
				    ++num_Ts;
				  }
				else if (type =="Tv")
				  {	//Tv
				    ++num_Tv;
				  }
			      }
		}
	    }
	  //make sure we use the right denominator
		
	  sites_compared = (ungapped_sites < seqlen) ? ungapped_sites : seqlen;
	  //P and Q are the proportions of Ts and Tv changes observed
	  P = double (num_Ts) / double (sites_compared);
	  Q = double (num_Tv) / double (sites_compared);
		
	  //Kimura's formula
	  if (fabs(1.0 - 2.0 * P - Q) > DBL_EPSILON)
	    {
	      divergence = -1.0 * 0.5 * log ((1.0 - 2.0 * P - Q)
	       *  pow ((1 - 2.0 * Q), 0.5));
	    }
	  else
	    {
	      divergence = 0.;
	    }
	  //a correction for extremely low observed values
	  if (divergence <= 0.0-DBL_EPSILON)
	    divergence = 0.0;
    }
              
      
      
    public:

       Kimura80 (const pair<string,string> *seqa, const pair<string,string> *seqb):
            seqlen (seqa->second.length())
       {
            if (seqa->second.length () != seqb->second.length ()){
                Rcpp::Rcout << "Kimura80(): constructor called with two sequence";
		Rcpp::Rcout << " objects of unequal length"<<endl;
                Rcpp::Rcout << seqa->second.length () << " vs. ";
		Rcpp::Rcout << seqb->second.length () <<endl;
            }
            num_Ts = 0;
            num_Tv = 0;
            divergence = 0.0;
            P = 0.0;
            Q = 0.0;
            sites_compared = 0;
            Compute (seqa,seqb);
       }
  
       double K (void)
        /*
          \return the distance between the two sequences.
          \note 999.0 is returned as a warning value.
          This can be necessary if sites are saturated,
          which implies that divergence cannot be calculated
        */
        {
            if (!isfinite(divergence))
              return (999.0);
            return (divergence);
        }

        size_t sites (void)
        /*
          \return the number of sites compared, excluding gaps, missing data, etc.
        */
        {
              return (sites_compared);
        }
};
    
