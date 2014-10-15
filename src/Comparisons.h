#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <assert.h>
#include <string>

#include "Sequence.h"

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


    
//    Ask if two strings are different.  While this can normally be done by asking
//    if (seq1 != seq2) {}, missing data poses a problem here.  If skip-missing == 1,
//    missing data (the 'N' character for nucleotide data, 'X' for amino acid)
//    are not used to determine if the sequences are different.  If nucleic_acid ==1,
//    nucleotide data are assumed, if nucleic_acid==0, protein data are assumed.  
//    \note case-insensitive
//    \return true if the seqs are different, false otherwise.  If the two sequences
//    are of different length, true is returned.
//    
// [[Rcpp::export]]
bool Different (const string & seq1,
        	    const string & seq2,
		    bool skip_missing = false,
		    bool nucleic_acid = false)  {
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
	      if( (_ch1 != MISSING && _ch2 != MISSING) && 
		  _ch1 != _ch2 )
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
unsigned NumDiffs (const string & seq1,
        	       const string & seq2,
		       bool skip_missing = false,
		       bool nucleic_acid = false)
    /*
      \return the number of differences between two strings.  Can skip missing
      data in the same fashion as Comparisons::Different.  If one sequence is shorter
      than the other, the number of positions compared is the length of the shorter 
      sequence.
    */
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
	      if( (_ch1 != MISSING && 
		   _ch2 != MISSING) && 
		  (_ch1 != _ch2 ) )
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
    
    


    
    //    Takes two ints, assumed to be integer representations of nucleotides. 
//    The way to ensure that the int represents a nucleotide in a valid way is
//    to use Sequence::Nucleotides.
//    The return value is determined
//    by a call to Comparisons::TsTv(int i, int j), where the ints are defined
//    in turn by Sequence::Nucleotides

    
// [[Rcpp::export]]
string TsTv (int i, int j) {
      //assert(i<=3 && j <= 3);
      if( ! (i<=3 && j <= 3)){ 
		cout << "Wrong Nucleotide "<<"i="<<i<<" j="<<j<<endl;
		cout << "Trying to fix ..."; 
		// whereever it comes from, there is still a call of TsTv without using nucToInt. 
		// We try to cath those here
		i = nucToInt( (char) i);
		j = nucToInt( (char) j);
		if( ! (i<=3 && j <= 3)){
			cout << "Still not right."<<endl;
			return "Unknown";
		}
		cout << endl;
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

      cout << "Neither Transition, not Transversion" <<endl;
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

