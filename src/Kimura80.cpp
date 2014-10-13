#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include <Comparisons.cpp>
// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar



class Kimura80
    {
    private:
      unsigned num_Ts, num_Tv;
      size_t seqlen;        	//total sequence length
      size_t sites_compared;	//number of ungapped sites in the data
      
      double divergence, P, Q;
      
      
      void Compute (const pair<string,string> *seq1, const pair<string,string> *seq2)
          {
            unsigned i;
        
            unsigned ungapped_sites = 0;
            for (i = 0; i < seqlen; ++i)        //iterate over the sequence
              {
                string type = "";
                //if( (strcmp(seq1->second[i],"-") != 0) && (strcmp(seq2->second[i],"-") != 0))
                if (NotAGap(seq1->second[i]) && NotAGap(seq2->second[i]))
                  {
                    ++ungapped_sites;
                    if (std::toupper(seq1->second[i]) != 
        		std::toupper(seq2->second[i]))	//if the sites differ at that position
                      {
                        type = TsTv (seq1->second[i], seq2->second[i]);	//check if difference is Ts or Tv
                      
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
              /*
                \param seqa an object of type Sequence::Seq
                \param seqb an object of type Sequence::Seq
                \exception Sequence::SeqException if sequences are of different lengths
              */
          {
            if (seqa->second.length () != seqb->second.length ()){
             // throw SeqException ("Sequence::Kimura80::Kimura80(): constructor called with two sequence objects of unequal length");
            cout << "Sequence::Kimura80::Kimura80(): constructor called with two sequence objects of unequal length"<<endl;
            cout << seqa->second.length () << " vs. " << seqb->second.length () <<endl;
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
    
