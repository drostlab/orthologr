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

#include <cassert>
#include <map>
#include <cctype>
#include <cstring>
#include <algorithm>
#include <functional>

#include "GranthamWeights.h"

class RedundancyCom95{
        
   //counting routines, return values from private matrices of same name
    
   
    double firstNon[4][4][4];
    double first2S[4][4][4];
    double first2V[4][4][4];
    double thirdFour[4][4][4];
    double thirdNon[4][4][4];
    double third2S[4][4][4];
    double third2V[4][4][4];
    
    double l0_vals[4][4][4];
    double l2S_vals[4][4][4];
    double l2V_vals[4][4][4];
    double l4_vals[4][4][4];
        

  void FillFirstPositionCounts ()
  {
    int i, j, k, l;
    string codon,mutation;
    //numbers of changes of various types, N = nonsynon, S = synon
    int numN, numTsS, numTvS;   
    int codonState;
    string type;		//type refers to Ts or Tv
    int numPossChanges, numPossTs, numPossTv;
    double FirstN, First2S, First2V;
    string codonTrans;
    string mutationTrans;

    bool S, N;
    //make every codon
    codon.resize(3);
    mutation.resize(3);
    for (i = 0; i <= 3; i++)
      {
        for (j = 0; j <= 3; j++)
          {
            for (k = 0; k <= 3; k++)
              {
                numN = numTsS = numTvS = 0;	//initialize
                numPossChanges = numPossTs = numPossTv = 0;
                S = N = 0;
                codon[0] = intToNuc(i);
                codon[1] = intToNuc(j);
                codon[2] = intToNuc(k);
  

                //now make all possible mutations at first position
                codonState = i;
                for (l = 0; l <= 3; l++)
                  {
                    if (l != codonState)
                      {	//only include mutations
                        mutation[0] = intToNuc(l);
                        mutation[1] = codon[1];
                        mutation[2] = codon[2];	
			//the mutant codon is now defined
                        type = TsTv( codonState, l);
			//find out if change is Ts or Tv

                        codonTrans = TranslateCodon (codon);
                        mutationTrans = TranslateCodon (mutation);

                        //if neither codon is a stop
                        if ((codonTrans)[0] != '*'
                            && (mutationTrans)[0] != '*')
                          {
                            if(codonTrans == mutationTrans)
                              {	//synonymous change
                                S = 1;
                                N = 0;
                              }
                            else
                              {	//nonsynonymous change
                                S = 0;
                                N = 1;
                              }
                            if (type == "Ts")
                              ++numPossTs;
                            else if (type == "Tv")
                              ++numPossTv;
                            if (N == 1)	//if change is nonsynonymous
                              ++numN;
                            else if (S == 1 && type == "Ts")
                              //if change is synonymous transition
                              ++numTsS;
                            else if (S == 1 && type == "Tv")
                              //if change is synonymous transversion
                              ++numTvS;
                          }
                      }
                  }
                //process here
                if ((numPossTs + numPossTv) == 0)
                  {	//stop codons
                    First2S = 0;
                    First2V = 0;
                    FirstN = 0;
                  }
                else if ((numTsS + numTvS) == 0)	//non-degenerate sites
                  {
                    First2S = 0;
                    First2V = 0;
                    FirstN = 1.0;
                  }
                else if (numPossTs != 0 && numPossTv != 0)
                  {	//special cases
                    if (double (numTsS) /
                        double (numPossTs) != 1.0
                        && double (numTvS) /
                        double (numPossTv) != 1.0)
                      {	//fractional redundancy
                        First2S =
                          double (numTsS) /
                          double (numPossTs +
                                  numPossTv);
                        First2V =
                          double (numTvS) /
                          double (numPossTs +
                                  numPossTv);
                        FirstN = 1.0 - First2S -
			  First2V;
                      }
                    else
                      {	//odd types of degeneracy
                        First2S =
                          double (numTsS) /
                          double (numPossTs);
                        First2V =
                          double (numTvS) /
                          double (numPossTv);
                        FirstN = 1.0 - First2S -
			  First2V;
                      }
                  }
                else
                  {
                    if (numPossTs > 0)
                      First2S =
                        double (numTsS) /
                        double (numPossTs + numPossTv);
                    else
                      First2S = 0.0;

                    if (numPossTv > 0)
                      First2V =
                        double (numTvS) /
                        double (numPossTv + numPossTs);
                    else
                      First2V = 0.0;
                    FirstN = 1.0 - First2S - First2V;
                  }

                //fill array
		firstNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=FirstN;
		first2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=First2S;
                first2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=First2V;

          }
        }
      }
  }


void FillThirdPositionCounts ()
  {
    int i, j, k, l;
    std::string codon, mutation;        //the two codons
    //numbers of changes of various types, N = nonsynon, S = synon
    int numN, numTsS, numTvS;	
    int codonState;
    string type;		//type refers to Ts or Tv
    int numPossChanges, numPossTs, numPossTv;
    double ThirdN, Third2S, Third2V, Third4;
    string codonTrans;
    string mutationTrans;
    bool S, N;
    //make every codon
    codon.resize(3);
    mutation.resize(3);
    for (i =0; i <= 3; i++)
      {
        for (j =0; j <= 3; j++)
          {
            for (k = 0; k <= 3; k++)
              {
                numN = numTsS = numTvS = 0;	//initialize
                numPossChanges = numPossTs = numPossTv = 0;
                S = N = 0;
                codon[0] = intToNuc(i);
                codon[1] = intToNuc(j);
                codon[2] = intToNuc(k);
                //now make all possible mutations at first position
                codonState = k;

                for (l = 0; l <= 3; ++l)
                  {
                    if (l != codonState)
                      {	//only include mutations
                        mutation[0] = codon[0];
                        mutation[1] = codon[1];
                        mutation[2] = intToNuc(l);	
			//the mutant codon is now defined
                        type = TsTv (codonState, l);	
			//find out if change is Ts or Tv

                        codonTrans = TranslateCodon (codon);
                        mutationTrans = TranslateCodon (mutation);

                        if ((codonTrans)[0] != '*'
                            && (mutationTrans)[0] != '*')
                          {
                            if(codonTrans==mutationTrans)
                              {
                                S = 1;
                                N = 0;
                              }
			    else
                              {
                                S = 0;
                                N = 1;
                              }
                            if (type == "Ts")
                              ++numPossTs;
                            else if (type == "Tv")
                              ++numPossTv;
                            if (N)
                              ++numN;
                            else if (S==1 && type == "Ts")
                              ++numTsS;
                            else if (S==1 && type == "Tv")
                              ++numTvS;
                          }
                      }
                  }
                //process here
                ThirdN = Third2S = Third2V = Third4 = 0.0;
                if (numPossTs + numPossTv == 0)
                  {	//stop codons
                    ThirdN = 0.0;
                    Third2S = 0.0;
                    Third2V = 0.0;
                    Third4 = 0.0;
                  }
                else if (numTsS + numTvS == 0)
                  {	//nondegenerate
                    ThirdN = 1.0;
                    Third2S = 0.0;
                    Third2V = 0.0;
                    Third4 = 0.0;
                  }
                else if (numTsS + numTvS == 3)
                  {	//4fold degenerate
                    ThirdN = 0.0;
                    Third2S = 0.0;
                    Third2V = 0.0;
                    Third4 = 1.0;
                  }
                else if (numTsS == 0 || numTvS == 0)
                  {	//if there is non-degeneracy for one type
                    if (numPossTs > 0)
                      Third2S =
                        double (numTsS) /
                        double (numPossTs);
                    else
                      Third2S = 0.0;
                    if (numPossTv > 0)
                      Third2V =
                        double (numTvS) /
                        double (numPossTv);
                    else
                      Third2V = 0.0;
                    ThirdN = 1.0 - Third2S - Third2V;
                    Third4 = 0.0;
                  }
                else if (numTsS > 0 && numTvS > 0
                         && (double (numTsS) /
                             double (numPossTs) != 1.0
                             || double (numTvS) /
                             double (numPossTv) != 1.0))
                  {	//fractional degeneracy
                    Third2S = 1.0 / 3.0;
                    Third2V = 1.0 / 3.0;
                    ThirdN = 1.0 / 3.0;
                    Third4 = 0.0;
                  }
		thirdNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=ThirdN;
                third2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=Third2S;
                third2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=Third2V;
		thirdFour[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])]=Third4;
              }
          }
      }
  }
  
  //the L values are just the sums of the site degeneracy values for each codon
  void FillLValues ()
  {
    int i, j, k;
    std::string trans;

    std::string codon(3,'A');
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
        for (k = 0; k <= 3; k++)
          {
            codon[0] = intToNuc(i);
            codon[1] = intToNuc(j);
            codon[2] = intToNuc(k);

            trans =TranslateCodon (codon);
            if(strcmp(trans.c_str(),"*")!=0)
           
            {        //if not a stop codon
		l0_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 1. 
		  + firstNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] 
		  + thirdNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];

		l2S_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 
		  first2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] + 
		  third2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];

		l2V_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] =
		  first2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] + 
		  third2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];

		l4_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 
		thirdFour[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];

                }
            else
              {
		l0_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 0.0;
		l2S_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 0.0;
		l2V_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 0.0;
		l4_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])] = 0.0;
              }
          }
  }
  

  public: 

  RedundancyCom95 ()
  {
    FillFirstPositionCounts ();
    FillThirdPositionCounts ();
    FillLValues ();    
  }
     
//the functions below are klugdy, but there is no real loss of efficieny,
//so it is left alone for now

  double FirstNon (const std::string & codon) const 
  /*
    \return number of times the first codon position is non-degenerate
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) {
        return firstNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    }
    Rcpp::Rcout << "Error in FirstNon - no codon"<< endl;
    return -1;
  }

  double First2S (const std::string & codon) const 
  /*
    \return number of times the first codon position is synonymous via a transition
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return first2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in First2S - no codon"<< endl;
    return -1;
  }

  double First2V (const std::string & codon) const 
  /*
    \return number of times the first codon position is synonymous via a transversion
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return first2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in First2V - no codon"<< endl;
    return -1;
  }

  double ThirdNon (const std::string & codon) const 
  /*
    \return number of times the third position is non-degenerate
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return thirdNon[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in ThirdNon - no codon"<< endl;
    return -1;
  }

  double ThirdFour (const std::string & codon) const 
  /*
    \return number of times the third position is fourfold-degenerate
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return thirdFour[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in ThirdFour - no codon"<< endl;
    return -1;
  }

  double Third2S (const std::string & codon) const 
  /*
    \return number of times the third position is synonymous via a transition
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return third2S[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in Third2s - no codon"<< endl;
    return -1;
  }

  double Third2V (const std::string & codon) const 
  /*
    \return number of times the third position is synonymous via a transversion
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
  */
  {
    if(codonPrecondition(codon)) 
        return third2V[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in Third2V - no codon"<< endl;
    return -1;
  }

  double L0_vals (const std::string & codon) const 
  /*
    \return the number of non-degenerate positions in \a codon
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
    \note the return value = 1.0+FirstNon(\a codon )+ThirdNon(\a codon)
  */
  {
    if(codonPrecondition(codon)) 
        return l0_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in L0_vals - no codon"<< endl;
    return -1;
  }

  double L2S_vals (const std::string & codon) const 
  /*
    \return the number of transitional silent sites in \a codon
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
    \note the return value = First2S(\a codon )+Third2S(\a codon )
  */
  {
    if(codonPrecondition(codon)) 
        return l2S_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in L2S_vals - no codon"<< endl;
    return -1;
  }

  double L2V_vals (const std::string & codon) const 
  /*
    \return the number of transversional silent sites in \a codon
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
    \note the return value = First2V(\a codon )+Third2V(\a codon )
  */
  {
    if(codonPrecondition(codon)) 
        return l2V_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    Rcpp::Rcout << "Error in L2V_vals - no codon"<< endl;
    return -1;
  }

  double L4_vals (const std::string & codon) const 
  /*
    \return the number of fourfold silent sites in \a codon
    \pre codon is of length 3, is all uppercase, and only contains the characters {A,G,C,T}
    \note the return value = ThirdFour(\a codon )
  */
  {
    if(codonPrecondition(codon)){
        return l4_vals[nucToInt(codon[0])][nucToInt(codon[1])][nucToInt(codon[2])];
    }
    Rcpp::Rcout << "Error in L4_vals - no codon"<< endl;
    return -1;
  }
  
};
