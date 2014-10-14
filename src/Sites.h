#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//#include <RedundancyCom95.cpp>
//#include <ThreeSubs.cpp>
#include "ThreeSubs.h"

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

class Sites {
        

      int seqlen;
      int maxhits;

      double _L0;
      double _L2S;
      double _L2V;
      double _L4;
      
      
                void siteinc (const RedundancyCom95 * sitesObj,
                          const std::string & codon1,
                          const std::string & codon2)
          {
            //skip ambiguous
        
            if ( !ambigousNucleotides(codon1) && !ambigousNucleotides(codon2) )
//                    std::find_if(codon1.begin(),codon1.end(),ambiguousNucleotide()) == codon1.end() &&
//                 std::find_if(codon2.begin(),codon2.end(),ambiguousNucleotide()) == codon2.end() )
              {
                _L0 += (sitesObj->L0_vals(codon1) + sitesObj->L0_vals(codon2))/2.0;
           //     cout << sitesObj->L0_vals(codon1) << "+" <<sitesObj->L0_vals(codon2) << "/2 =" << _L0<<endl;
                _L2S += (sitesObj->L2S_vals(codon1) + sitesObj->L2S_vals(codon2))/2.0;
                _L2V += (sitesObj->L2V_vals(codon1) + sitesObj->L2V_vals(codon2))/2.0;
                _L4 += (sitesObj->L4_vals(codon1) + sitesObj->L4_vals(codon2))/2.0;
              }
          } 

        
        void count_sites (const pair<string,string> * sequence1,
                      const pair<string,string> * sequence2,
                      const RedundancyCom95 * sitesObj)
          {        			//now, it is correct
            size_t i = 0, j = 0;
            int nc;
            std::string codon1,codon2;
            codon1.resize(3);
            codon2.resize(3);
            for (i = 0; i <= seqlen - 3; i += 3)
              {
                for (j = 0; j <= 2; ++j)
                  {
                    codon1[j] = char(std::toupper((sequence1->second)[i + j]));
                    codon2[j] = char(std::toupper((sequence2->second)[i + j]));
                  }
                nc = NumDiffs (codon1, codon2);
//        cout << "Number of diffs between "<< codon1 << " and " << codon2 << " are " <<nc<<endl;
                if (nc == 0)	//still need to count if there are 0 changes
                  siteinc (sitesObj,codon1, codon2);
                else if (maxhits <= 3 && nc == 1)
                  siteinc (sitesObj,codon1, codon2);
                else if (maxhits == 2 && nc <= 2)
                  siteinc (sitesObj,codon1, codon2);
                else if (maxhits == 3 && nc <= 3)
                  siteinc (sitesObj,codon1, codon2);
              }
          }


          
                public:
        Sites (const RedundancyCom95 * sitesObj,const pair<string,string> * seq1,
                                  const pair<string,string> * seq2, 
                                  const int max  ):
              _L0(0.),_L2S(0.),_L2V(0.0),_L4(0.0)
              /*
                \param sitesObj an initialized object of type Sequence::RedundancyCom95
                \param seq1 a Sequence::Seq
                \param seq2 a Sequence::Seq
                \param max max number of substitutions per codon to analyze
                \param code see Sequence::GeneticCodes for valid values
                \note sequences must be of same length, this is checked by assert()
                \note sequence lengths must be multiples of 3, this is checked by assert()
               */
          {
            assert(seq1->second.length() == seq2->second.length());
            assert( (fabs(seq1->second.length()%3-0.)<=DBL_EPSILON) &&
                    (fabs(seq1->second.length()%3-0.)<=DBL_EPSILON ) );
            seqlen = seq1->second.length();
            maxhits = max;
          //  genetic_code = code;
            count_sites (seq1,seq2,sitesObj);
          }
          
                double L0(void) const
      /*!
        \return alignment length in terms of non-degenerate sites
      */
        {
          return _L0;
        }
      double L2S(void) const
      /*!
        \return alignment length in terms of transitional-degenerate sites
      */
        {
          return _L2S;
        }
      double L2V(void) const
      /*!
        \return alignment length in terms of transversional-degenerate sites
      */
        {
          return _L2V;
        }
      double L4(void) const
      /*!
        \return alignment length in terms of fourfold-degenerate sites
      */
        {
          return _L4;
        }
         
};




