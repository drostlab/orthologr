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

#include <algorithm>
#include <cctype>

// [[Rcpp::export]]
bool ambigousNucleotides(const std::string & codon){
        for(int i=0;i<3;i++){
                if(! ((codon)[i] == 'A' 
                || (codon)[i] == 'C' 
                || (codon)[i] == 'G' 
                || (codon)[i] == 'T' ))
                { return true; }
        }
        return false;
}

  
// [[Rcpp::export]]
bool codonPrecondition(const std::string & codon)
{
    if ( codon.length() != 3 || !ambigousNucleotides(codon) )
         return true;
    return false;
}
 

// [[Rcpp::export]]
char intToNuc(int i)
{
          switch (i){
                  case 0: return 'A';break;
                  case 1: return 'T';break;
                  case 2: return 'G';break;
                  case 3: return 'C';break;
          }
          return 'X';
}
  

// [[Rcpp::export]]
int nucToInt(char c)
{
//the map is case-insensitive...
          switch (c){
                  case 'A': return 0;break;
                  case 'T': return 1;break;
                  case 'G': return 2;break;
                  case 'C': return 3;break;
                  case 'a': return 0;break;
                  case 't': return 1;break;
                  case 'g': return 2;break;
                  case 'c': return 3;break;
		  //default: return -1;
          }
          return -1;
}
 

// [[Rcpp::export]]
char Universal(std::string codon)
  {
    //handle gaps A codon with a single or 2
    //gap characters returns an X, since it's ambiguous,
    //and a codon that is all gaps (---)
    //returns a -
        char gapchar='-';
      int ngaps = (count(&codon[0],&codon[0]+3,gapchar));
        
     if (ngaps == 3)
      {
        return '-';
      }
    else if (ngaps > 0 && ngaps < 3)
      {
	return 'X';
      }

    switch (codon[0])
      {
        //first position is A
      case 'A':
        switch (codon[1])
          {
          case 'T':
            switch (codon[2])
              {
              case 'G':
                return 'M';
                break;
              default:
                return 'I';
                break;
              }
            break;
          case 'C':
            return 'T';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'T':
                return 'N';
                break;
              case 'C':
                return 'N';
                break;
              case 'A':
                return 'K';
                break;
              case 'G':
                return 'K';
                break;
              }
          case 'G':
            switch (codon[2])
              {
              case 'A':
                return 'R';
                break;
              case 'G':
                return 'R';
                break;
              default:
                return 'S';
                break;
              }
            break;
          }
        break;
        //first Position is T
      case 'T':
        switch (codon[1])
          {
          case 'T':
            switch (codon[2])
              {
              case 'T':
                return 'F';
                break;
              case 'C':
                return 'F';
                break;
              case 'A':
                return 'L';
                break;
              case 'G':
                return 'L';
                break;
              }
            break;
          case 'C':
            return 'S';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'T':
                return 'Y';
                break;
              case 'C':
                return 'Y';
                break;
              default:
                return '*';
                break;
              }
            break;
          case 'G':
            switch (codon[2])
              {
              case 'A':
                return '*';
                break;
              case 'G':
                return 'W';
                break;
              default:
                return 'C';
                break;
              }
            break;
          }

        break;
        //first is G
      case 'G':
        switch (codon[1])
          {
          case 'T':
            return 'V';
            break;
          case 'C':
            return 'A';
            break;
          case 'G':
            return 'G';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'A':
                return 'E';
                break;
              case 'G':
                return 'E';
                break;
              default:
                return 'D';
                break;
              }
            break;
          }
        break;
        //first is C
      case 'C':
        switch (codon[1])
          {
          case 'T':
            return 'L';
            break;
          case 'C':
            return 'P';
            break;
          case 'G':
            return 'R';
            break;
          case 'A':
            switch (codon[2])
              {
              case 'A':
                return 'Q';
                break;
              case 'G':
                return 'Q';
                break;
              default:
                return 'H';
                break;
              }
            break;
          }
        break;
      default:
        return 'X';
        break;
      }
    return 'X';
 }

 
// [[Rcpp::export]]
std::string TranslateCodon(std::string codon)
  {
    //if the range is less than 3 in length (1 codon), return an empty string
    if (codon.length() != 3) 
        return std::string();

    std::string translation;

        codon[0] = char(std::toupper(codon[0]));
        codon[1] = char(std::toupper(codon[1]));
        codon[2] = char(std::toupper(codon[2]));

        translation += Universal (codon);

    return translation;
  }
