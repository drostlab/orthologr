#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include <algorithm>
#include <cctype>

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
bool ambigousNucleotides(const std::string & codon){
        for(int i=0;i<3;i++){
//                cout << codon[i] << endl;
                if(! ((codon)[i] == 'A' 
                || (codon)[i] == 'C' 
                || (codon)[i] == 'G' 
                || (codon)[i] == 'T' ))
                { return true; }
        }
        return false;
}

//struct ambiguousNucleotide : public std::unary_function<char,bool>
//                               /*! \struct ambiguousNucleotide Sequence/SeqProperties.hpp
//				*/
//  {
//    inline bool operator()(const char & c) const
//    /*!
//      \return true if c is not A,G,C, or T, false otherwise
//      \note Case-insensitive
//    */
//    {
//      const char ch = char(std::toupper(c));
//      return (ch != 'A' &&
//	      ch != 'G' &&
//	      ch != 'T' &&
//	      ch != 'C' );
//    }
//  };
  
  
  /*
    Checks to see if this class can handle the codon.
    A valid codon is:
    1.) of length 3, 
    2.) contains only characters in the set {A,G,C,T} (case sensitive)
  */  
// [[Rcpp::export]]
bool codonPrecondition(const std::string & codon)
  {
  //  cout << "Precondition length = " << codon.length() << endl;
     
    
    if ( codon.length() != 3 || !ambigousNucleotides(codon) )
     //    std::find_if(codon.begin(),codon.end(),ambiguousNucleotide()) != codon.end() )
         return true;
        return false;
  }
 

  
// [[Rcpp::export]]
char intToNuc(int i){
          switch (i){
                  case 0: return 'A';break;
                  case 1: return 'T';break;
                  case 2: return 'G';break;
                  case 3: return 'C';break;
          }
  }
  
    //the map is case-insensitive...
// [[Rcpp::export]]
int nucToInt(char c){
          switch (c){
                  case 'A': return 0;break;
                  case 'T': return 1;break;
                  case 'G': return 2;break;
                  case 'C': return 3;break;
                  case 'a': return 0;break;
                  case 't': return 1;break;
                  case 'g': return 2;break;
                  case 'c': return 3;break;
          }
  }
 

// [[Rcpp::export]]
char Universal(string codon)
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
std::string TranslateCodon(string codon)
  {
//             cout << "Translate "<<codon;

    if (codon.length() != 3) //if the range is less than 3 in length (1 codon), return an empty string
      return std::string();

    string translation;

        codon[0] = char(std::toupper(codon[0]));
        codon[1] = char(std::toupper(codon[1]));
        codon[2] = char(std::toupper(codon[2]));
//                cout << "  Universal of  "<<codon<<endl;

            translation += Universal (codon);

 
    return translation;
  }