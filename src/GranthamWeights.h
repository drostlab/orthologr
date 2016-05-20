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
#include "Grantham.h"
#include "Kimura80.h"

class GranthamWeights2
{
  public:  

  GranthamWeights2()
  /*
        \param genetic_code the code used to TranslateCodon the codons
  */
  {
  }
  
  void Intermediates2( string *intermediates, string codon1, string codon2)
  {  
    intermediates[0].resize(3);
    intermediates[1].resize(3);

    unsigned i,j;
    int pos[2];

    for (i = 0, j = 0; i <= 2; ++i)
      if (char(toupper(codon1[i])) != char(toupper(codon2[i])))
        pos[j++] = i;
    switch (pos[0])
      {
      case 0:
	intermediates[0][0] = char(toupper(codon2[0]));
	intermediates[0][1] = char(toupper(codon1[1]));
	intermediates[0][2] = char(toupper(codon1[2]));
	break;
      case 1:
	intermediates[0][0] = char(toupper(codon1[0]));
	intermediates[0][1] = char(toupper(codon2[1]));
	intermediates[0][2] = char(toupper(codon1[2]));
	break;
      case 2:
	intermediates[0][0] = char(toupper(codon1[0]));
	intermediates[0][1] = char(toupper(codon1[1]));
	intermediates[0][2] = char(toupper(codon2[2]));
	break;
      }

    switch(pos[1])
      {
      case 0:
	intermediates[1][0] = char(toupper(codon2[0]));
	intermediates[1][1] = char(toupper(codon1[1]));
	intermediates[1][2] = char(toupper(codon1[2]));
	break;
      case 1:
	intermediates[1][0] = char(toupper(codon1[0]));
	intermediates[1][1] = char(toupper(codon2[1]));
	intermediates[1][2] = char(toupper(codon1[2]));
	break;
      case 2:
	intermediates[1][0] = char(toupper(codon1[0]));
	intermediates[1][1] = char(toupper(codon1[1]));
	intermediates[1][2] = char(toupper(codon2[2]));
	break;
      }
 }
        
  void Calculate(double *__weights, const string &codon1, const string &codon2){
  /*
    Calculate actually calculates the weights for each branch
    \param codon1 a string of length 3 representing a sense codon
    \param codon2 a string of length 3 representing a sense codon
  */

      Grantham gdist;
      string intermediates[2];
      Intermediates2(intermediates, codon1,codon2);
      //assign weights to pathways
      //weights are assiged by the total length of each path, as
      //measured by Grantham distances
      double len_path_1 = 0.0, len_path_2 = 0.0;

      string t1 = TranslateCodon (codon1);
      string t2 = TranslateCodon (intermediates[0]);
      len_path_1 += gdist ((t1)[0], (t2)[0]);

      t1 = TranslateCodon (intermediates[0]);
      t2 = TranslateCodon (codon2);
      len_path_1 += gdist ((t1)[0], (t2)[0]);

      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[1]);
      len_path_2 += gdist ((t1)[0], (t2)[0]);

      t1 = TranslateCodon (intermediates[1]);
      t2 = TranslateCodon (codon2);
      len_path_2 += gdist ((t1)[0], (t2)[0]);

      //calculate the weights themselves
      //double w_path1 = 0., w_path2 = 0.,
      double w_tot = 0.;
 
      __weights[0] = 0.;
      __weights[1]= 0.;

      if (fabs(len_path_1-0.) <= DBL_EPSILON && fabs(len_path_2-0.) <= DBL_EPSILON)
        {
          __weights[0] = __weights[1] = 0.5;
          w_tot = 1.;
        }
      else
        {
          __weights[0] = 1. - len_path_1 / (len_path_1 + len_path_2);
          __weights[1] = 1. - len_path_2 / (len_path_1 + len_path_2);
          w_tot = __weights[0] + __weights[1];
        }

      __weights[0] /= w_tot;
      __weights[1] /= w_tot;

    }
 
};

class GranthamWeights3{
                    
   public:

   GranthamWeights3()
   /*
        \param genetic_code the code used to TranslateCodon the codons
   */
   {
   }

  void Intermediates3( string *intermediates, string codon1,  string codon2)
  {
          
    for(int i = 0 ; i < 9 ;++i)
      intermediates[i].resize(3);

    intermediates[0][0] = char(toupper(codon2[0]));
    intermediates[0][1] = char(toupper(codon1[1]));
    intermediates[0][2] = char(toupper(codon1[2]));

    intermediates[1][0] = char(toupper(codon2[0]));
    intermediates[1][1] = char(toupper(codon2[1]));
    intermediates[1][2] = char(toupper(codon1[2]));

    intermediates[2][0] = char(toupper(codon2[0]));
    intermediates[2][1] = char(toupper(codon1[1]));
    intermediates[2][2] = char(toupper(codon2[2]));

    intermediates[3][0] = char(toupper(codon1[0]));
    intermediates[3][1] = char(toupper(codon2[1]));
    intermediates[3][2] = char(toupper(codon1[2]));

    intermediates[4][0] = char(toupper(codon2[0]));
    intermediates[4][1] = char(toupper(codon2[1]));
    intermediates[4][2] = char(toupper(codon1[2]));

    intermediates[5][0] = char(toupper(codon1[0]));
    intermediates[5][1] = char(toupper(codon2[1]));
    intermediates[5][2] = char(toupper(codon2[2]));

    intermediates[6][0] = char(toupper(codon1[0]));
    intermediates[6][1] = char(toupper(codon1[1]));
    intermediates[6][2] = char(toupper(codon2[2]));

    intermediates[7][0] = char(toupper(codon2[0]));
    intermediates[7][1] = char(toupper(codon1[1]));
    intermediates[7][2] = char(toupper(codon2[2]));

    intermediates[8][0] = char(toupper(codon1[0]));
    intermediates[8][1] = char(toupper(codon2[1]));
    intermediates[8][2] = char(toupper(codon2[2]));
  }

void Calculate(double *__weights, const string &codon1, const string &codon2) {
  /*
    Calculate actually calculates the weights for each branch
    \param codon1 a string of length 3 representing a sense codon
    \param codon2 a string of length 3 representing a sense codon
  */
    
      Grantham gdist;
      
      string intermediates[9];
      Intermediates3(intermediates,codon1,codon2);
      double len_path_1 = 0.0, len_path_2 = 0.0, len_path_3 = 0., 
        len_path_4 = 0., len_path_5 = 0., len_path_6 = 0.;
      double dist = 0.;

      string t1 = TranslateCodon (codon1);
      string t2 = TranslateCodon (intermediates[0]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;

      t1 = TranslateCodon (intermediates[0]);
      t2 = TranslateCodon (intermediates[1]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;

      t1 = TranslateCodon (intermediates[1]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_1 += dist;

      //path2

      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[0]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      t1 = TranslateCodon (intermediates[0]);
      t2 = TranslateCodon (intermediates[2]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      t1 = TranslateCodon (intermediates[2]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_2 += dist;

      //path 3
      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[3]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      t1 = TranslateCodon (intermediates[3]);
      t2 = TranslateCodon (intermediates[4]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      t1 = TranslateCodon (intermediates[4]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_3 += dist;

      //path4
      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[3]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      t1 = TranslateCodon (intermediates[3]);
      t2 = TranslateCodon (intermediates[5]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      t1 = TranslateCodon (intermediates[5]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_4 += dist;

      //path 5
      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[6]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      t1 = TranslateCodon (intermediates[6]);
      t2 = TranslateCodon (intermediates[7]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      t1 = TranslateCodon (intermediates[7]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_5 += dist;

      //path 6
      t1 = TranslateCodon (codon1);
      t2 = TranslateCodon (intermediates[6]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      t1 = TranslateCodon (intermediates[6]);
      t2 = TranslateCodon (intermediates[8]);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      t1 = TranslateCodon (intermediates[8]);
      t2 = TranslateCodon (codon2);
      dist = gdist ((t1)[0],(t2)[0]);
      len_path_6 += dist;

      __weights[0] = 1. / len_path_1;
      __weights[1] = 1. /len_path_2;
      __weights[2] = 1. /len_path_3;
      __weights[3] = 1. / len_path_4;
      __weights[4] = 1. /len_path_5;
      __weights[5] = 1./ len_path_6;

      //scale weights to sum to 1
      double w_tot = __weights[0] + __weights[1] + __weights[2] + __weights[3] + __weights[4] + __weights[5];
      __weights[0] /= w_tot;
      __weights[1] /= w_tot;
      __weights[2] /= w_tot;
      __weights[3] /= w_tot;
      __weights[4] /= w_tot;
      __weights[5] /= w_tot;
    }  

};
