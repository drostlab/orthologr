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

#include "SingleSub.h"

class TwoSubs{

  private:
      double p0, p2S, p2V, p4, q0, q2S, q2V,q4;
      double p0_b1, p2S_b1, p2V_b1, p4_b1, q0_b1, q2S_b1, q2V_b1, q4_b1;
      double p0_b2, p2S_b2, p2V_b2, p4_b2, q0_b2, q2S_b2, q2V_b2, q4_b2;
      double p0_b3, p2S_b3, p2V_b3, p4_b3, q0_b3, q2S_b3, q2V_b3, q4_b3;
      double p0_b4, p2S_b4, p2V_b4, p4_b4, q0_b4, q2S_b4, q2V_b4, q4_b4;
     
   void Calculate (const RedundancyCom95 * sitesObj, const std::string & codon1,
                   const std::string & int_1, const std::string & int_2,
                   const std::string & codon2, const double w_path1,
                   const double w_path2)
  /*
    count up mutations between the codons
  */
  {
    //count the changes along each branch
    //branches that go through stop codons are counted as
    //containing no changes on them.
    SingleSub Single;
    Single(sitesObj,codon1,int_1);
    p0_b1 = Single.P0();
    p2S_b1 =Single.P2S();
    p2V_b1 =Single.P2V();
    p4_b1 = Single.P4();
    q0_b1 = Single.Q0();
    q2S_b1 =Single.Q2S();
    q2V_b1 =Single.Q2V();
    q4_b1 = Single.Q4();

    Single(sitesObj, int_1, codon2);
    p0_b2 = Single.P0();
    p2S_b2 =Single.P2S();
    p2V_b2 =Single.P2V();
    p4_b2 = Single.P4();
    q0_b2 = Single.Q0();
    q2S_b2 =Single.Q2S();
    q2V_b2 =Single.Q2V();
    q4_b2 = Single.Q4();


     Single(sitesObj, codon1, int_2);
    p0_b3 = Single.P0();
    p2S_b3 =Single.P2S();
    p2V_b3 =Single.P2V();
    p4_b3 = Single.P4();
    q0_b3 = Single.Q0();
    q2S_b3 =Single.Q2S();
    q2V_b3 =Single.Q2V();
    q4_b3 = Single.Q4();

     Single(sitesObj, int_2, codon2);
    p0_b4 = Single.P0();
    p2S_b4 =Single.P2S();
    p2V_b4 =Single.P2V();
    p4_b4 = Single.P4();
    q0_b4 = Single.Q0();
    q2S_b4 =Single.Q2S();
    q2V_b4 =Single.Q2V();
    q4_b4 = Single.Q4();


    p0  = (p0_b1  + p0_b2)  * w_path1
          + (p0_b3  + p0_b4) * w_path2;
    p2S = (p2S_b1 + p2S_b2) * w_path1
          + (p2S_b3 + p2S_b4) * w_path2;
    p2V = (p2V_b1 + p2V_b2) * w_path1
          + (p2V_b3 + p2V_b4) * w_path2;
    p4  = (p4_b1  + p4_b2)  * w_path1
          + (p4_b3  + p4_b4) * w_path2;
    q0  = (q0_b1  + q0_b2)  * w_path1
          + (q0_b3  + q0_b4) * w_path2;
    q2S = (q2S_b1 + q2S_b2) * w_path1
          + (q2S_b3 + q2S_b4) * w_path2;
    q2V = (q2V_b1 + q2V_b2) * w_path1
          + (q2V_b3 + q2V_b4) * w_path2;
    q4  = (q4_b1  + q4_b2)  * w_path1
          + (q4_b3  + q4_b4) * w_path2;
  }



  public:
    
     void operator()(const RedundancyCom95 * sitesObj,
                     const std::string & codon1, const std::string & codon2,
                     GranthamWeights2 *weights2)
  /*
    \param sitesObj an initialized object of type RedundancyCom95
    \param code genetic code, see GeneticCodes for valid values
    \param codon1 string of length 3 representing nucleodites
    \param codon2 string of length 3 representing nucleodites
    \param weights2 a weighting scheme for the pathways
  */
  {
    p0= p2S= p2V= p4= q0= q2S= q2V=q4=0.;
    p0_b1= p2S_b1= p2V_b1= p4_b1= q0_b1= q2S_b1= q2V_b1= q4_b1=0.;
    p0_b2= p2S_b2= p2V_b2= p4_b2= q0_b2= q2S_b2= q2V_b2= q4_b2=0.;
    p0_b3= p2S_b3= p2V_b3= p4_b3= q0_b3= q2S_b3= q2V_b3= q4_b3=0.;
    p0_b4= p2S_b4= p2V_b4= p4_b4= q0_b4= q2S_b4= q2V_b4= q4_b4=0.;
    

    std::string intermediates[2];
    weights2->Intermediates2(intermediates,codon1,codon2);
    
    double weights[2];
    weights2->Calculate(weights,codon1,codon2);

    Calculate (sitesObj, codon1, intermediates[0], codon2, intermediates[1], weights[0], weights[1]); 
 }

  double P0 (void) const
  /*
    \return number of transitions at non-degenerate sites in the codon
  */
  {
    return p0;
  }

  double P2S (void) const
  /*
    \return number of transitions at transitional-degenerate sites in the codon
  */
  {
    return p2S;
  }

  double P2V (void) const
  /*
     \return number of transitions at transversional-degenerate sites in the codon
   */
  {
    return p2V;
  }

  double P4 (void) const
  /*
    \return number of transitions at fourfold-degenerate sites in the codon
  */
  {
    return p4;
  }

  double Q0 (void) const
  /*
    \return number of transversions at non-degenerate sites in the codon
  */
  {
    return q0;
  }

  double Q2S (void) const
  /*
    \return number of transversions at transitional-degenerate sites in the codon
  */
  {
    return q2S;
  }

  double Q2V (void) const
  /*
    \return number of transversions at transversional-degenerate sites in the codon
  */
  {
    return q2V;
  }

  double Q4 (void) const
  /*
    \return number of transversions at fourfold-degenerate sites in the codon
  */
  {
    return q4;
  }

};
