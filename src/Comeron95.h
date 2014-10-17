#include <Rcpp.h>
#include <assert.h>     /* assert */
#include <math.h> 
#include <string>
using namespace Rcpp;
using namespace std;

//#include <Comparisons.cpp>
//#include <Sequence.cpp>
//#include <RedundancyCom95.cpp>

#include "Sites.h"
//#include <GranthamWeights.cpp>
//#include <SingleSub.cpp>

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


class Comeron95{


private:

        double Qs, Bs, Qa, Ba, A2S, A4, As, A2V, A0, Aa;
        double p0, p2V, p2S, p4, q0, q2V, q2S, q4;
      Sites *sites;
     const  RedundancyCom95 *sitesObj;
        double dN;
        double dS;
        double N;
        double S;     // number of nonsynonymous and synonymous substitutions respectively
        int maxhits;
        


        
void diverge(const pair<string,string> *seq1, const pair<std::string,std::string> *seq2,
           GranthamWeights2 *_weights2 = NULL, GranthamWeights3 *_weights3 = NULL
                     ){
               /*
    go through every aligned, ungapped codon,
    and calculate divergence.  maintains a running sum of divergence
    statistics stored a private data to the class
  */
    size_t i, j, ndiff;
    size_t length = seq1->second.length();

    //the for loop iterates over codons (block of 3 sites)
    std::string codon1, codon2;
    codon1.resize (3);
    codon2.resize (3);

    for (i = 0; i <= (length - 3); i += 3)
      {
        for (j = 0; j <= 2; ++j)
          {
            //assign the next codon from each sequence
            codon1[j] = char(std::toupper(seq1->second[i + j]));
            codon2[j] = char(std::toupper(seq2->second[i + j]));
          }

      if ( !ambigousNucleotides(codon1) && !ambigousNucleotides(codon2) )
//                std::find_if(codon1.begin(),codon1.end(),ambiguousNucleotide()) == codon1.end() &&
//              std::find_if(codon2.begin(),codon2.end(),ambiguousNucleotide()) == codon2.end() )
          {
            //find out if codons are different
            if (Different(codon1, codon2))
              {
                //find out how many changes there are between the codons
                ndiff = NumDiffs (codon1, codon2);
 
                if (ndiff == 1)	//if there is just one difference, use the rules for
                  //codons that only differ at 1 site
                  {
                    SingleSub Single;
                    Single(sitesObj, codon1, codon2);
                    p0 += Single.P0 ();
                    p2S += Single.P2S ();
                    p2V += Single.P2V ();
                    p4 += Single.P4 ();
                    q0 += Single.Q0 ();
                    q2S += Single.Q2S ();
                    q2V += Single.Q2V ();
                    q4 += Single.Q4 ();
                  }

                if (maxhits > 1)
                  {	//if codons with >1 difference are allowed
                    //iterate over codons, as above

// redundant
//                if (Different(codon1, codon2))
//                      {
//                              ndiff = NumDiffs(codon1,codon2);

                        //cout << "CHECK:NDIFF="<<ndiff<<endl;
                        if (ndiff == 2
                            && maxhits >= 2)
                          {	//if codons differ at 2 sites,
                            //use rules of class TwoSubstitutions
                            TwoSubs Double;
                            Double(sitesObj, codon1, codon2, _weights2);
                            p0 += Double.P0();
                            p2S += Double.P2S ();
                            p2V += Double.P2V ();
                            p4 += Double.P4();
                            q0 += Double.Q0 ();
                            q2S += Double.Q2S ();
                            q2V += Double.Q2V ();
                            q4 += Double.Q4 ();
                          }
                        else if (ndiff == 3 && maxhits > 2)
                          {
                            ThreeSubs Triple;
                            Triple(sitesObj,codon1, codon2,_weights3);
                            p0 += Triple.P0 ();
                            p2S += Triple.P2S ();
                            p2V += Triple.P2V ();
                            p4 += Triple.P4 ();
                            q0 += Triple.Q0 ();
                            q2S += Triple.Q2S ();
                            q2V += Triple.Q2V ();
                            q4 += Triple.Q4 ();
                          }
                      //}
                  }
              }
          }
      }
      
}
        
void omega(const pair<string,string> *seq1, const pair<std::string,std::string> *seq2){
/*
    calculate values needed to obtain dN and dS.
    formulae are from Comeron '95 and use the identical notation
    the rest of the code is just calculating numbers from the data.
    At first glance, it loodS like this is a lot of code,
    but these calculations require careful checking, becuase when you
    have a lot of changes, or very few, you can take logs of negative
    numbers, or divide by zero, hence calls to isnan() and isinf()
  */
  
    double log1, log2;

    Qs = (q2V + q4) / (sites->L2V() + sites->L4());

    if (!isfinite (Qs))
      Qs = 0.0;

    Bs = (-0.5) * log (1.0 - (2.0 * Qs));

         if (isnan (Bs))
           {
             Kimura80 *K80 = new Kimura80 (seq1, seq2);
     	//if sites aren't saturated in general (i.e. Kimura's distance < 1.0)
     	//it is likely that Bs is nan due to too few changes, and thus Bs should equal 0.0
     	//otherwise it is due to too many changes.
     	//NOTE--this is an ad-hoc treatment of the analysis!!!
     	if (K80->K () < 1.0)
     	  Bs = 0.0;
     	//delete (K80);
           }

    Qa = (q0 + q2S) / (sites->L0() + sites->L2S());

    if (!isfinite (Qa))
      Qa = 0.0;

    Ba = (-0.5) * log (1.0 - (2.0 * Qa));
    if (!isfinite (Ba))// && !isinf(Ba))
      {
        //if sites aren't saturated in general (i.e. Kimura's distance < 1.0)
        //it is likely that Ba is nan due to too few changes, and thus Ba should equal 0.0
        //otherwise it is due to too many changes.
        //NOTE--this is an ad-hoc treatment of the analysis!!!
        Kimura80 *K80 = new Kimura80 (seq1, seq2);
        if (K80->K () < 1.0)
          Ba = 0.0;
        //delete (K80);
      }

    //calculate numbers of mutation per site type
    double P2S_site = p2S / sites->L2S();
    double P2V_site = p2V / sites->L2V();
    double P0_site = p0 / sites->L0();
    double Q0_site = q0 / sites->L0();
    double P4_site = p4 / sites->L4();
    double Q4_site = q4 / sites->L4();

    log1 = log (1.0 - (2.0 * P2S_site) - Qa);
    log2 = log (1.0 - (2.0 * Qa));

    if (!isfinite (log1))	//set value to 0 if nan
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A2S = (-0.5) * log1 + (0.25) * log2;

    log1 = log (1.0 - (2.0 * P4_site) - Q4_site);
    log2 = log (1.0 - (2.0 * Q4_site));

    if (!isfinite (log1))
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A4 = (-0.5) * log1 + (0.25) * log2;

    As = (sites->L2S() * A2S + sites->L4() * A4) / (sites->L2S() +
         sites->L4());

    log1 = log (1.0 - (2.0 * P2V_site) - Qs);
    log2 = log (1.0 - (2.0 * Qs));

    if (!isfinite (log1))
      log1 = 0.0;
    if (!isfinite (log2))
      log2 = 0.0;

    A2V = (-0.5) * log1 + (0.25) * log2;

    log1 = log (1.0 - (2.0 * P0_site) - Q0_site);
    log2 = log (1.0 - (2.0 * Q0_site));
    if (!isfinite (log1))
      log1 = 0.;
    if (!isfinite (log2))
      log2 = 0.0;

    A0 = (-0.5) * log1 + (0.25) * log2;

    Aa = (sites->L2V() * A2V + sites->L0() * A0) / (sites->L2V() +
         sites->L0());

    if (As <= 0.0)
      As = 0.0;
    if (Bs <= 0.0)
      Bs = 0.0;
    if (Aa <= 0.0)
      Aa = 0.0;
    if (Ba <= 0.0)
      Ba = 0.0;

    dS = As + Bs;
    dN = Aa + Ba;

    if (!isfinite (dS))
      dS = 999;
    if (!isfinite (dN))
      dN = 999;
  }


public:
        Comeron95(const pair<string,string> *seq1, const pair<string,string> *seq2, int max, 
                 const RedundancyCom95 *genetic_code_redundancy){
                
                // max in 1-3
                assert(max>=1 && max <=3);
                maxhits = max;
                
                // proof if seqs have same length
                if( seq1->second.length() != seq2->second.length() ){
                        cout << "Sequences of unequal length.\n";
                }
                
                // proof if seqs are multiple of 3
                if(! ( (seq1->second.length()%3 == 0) && (seq2->second.length()%3 == 0)) ){
                        cout << "Sequence lengths are not multiples of 3. ";
                        cout << seq1->second.length() << " and " << seq2->second.length() <<endl;
                }
                
                // set weighting schemes
                
//            if(_weights2==NULL)
//              {
//                __2wasNULL=true;
               GranthamWeights2 * weights2 = new GranthamWeights2();
//              }
        
//            if( _weights3 == NULL)
//              {
//                __3wasNULL=true;
                GranthamWeights3 * weights3 = new GranthamWeights3();
//              }
        
                
                
                // prepare calculation building a Redundancy95 Object storing the number of 
                // mutations that can occur and how to be count
                
                 //If no weighting schemes are passed to the constructor,
    //then we default to using Grantham's distances

    if (genetic_code_redundancy == NULL)
      {
        sitesObj = new RedundancyCom95();
        //__red_was_NULL = true;
      }
    else
      sitesObj = genetic_code_redundancy;                     
     

    
                // new Sites Object counting L0-L4 of the sequence
      sites = new Sites(sitesObj, seq1, seq2, maxhits);

                q0=q2S=q2V=q4=p0=p2S=p2V=p4=0.0;
                
                diverge(seq1, seq2
                        //weights2 &3
                        );
                        
                omega( seq1, seq2);
}

        double dn (void) const
           /*
            \return the nonsynonymous distance
            \note 999.0 is returned if dN cannot be calculated
          */
          {
                return dN;
          }
        double ds (void) const
           /*
            \return the nonsynonymous distance
            \note 999.0 is returned if dN cannot be calculated
          */
          {
                return dS;
          }
        double ratio (void) const
           /*
            \return the nonsynonymous distance
            \note 999.0 is returned if dN cannot be calculated
          */
          {
                  // if dN==999 || dS==999
                return dN/dS;
          }
        
        
        
        
          double P0 (void) const
          /*
            \return number of transitions at nondegenerate sites
          */
          {
            return p0;
          }
          double P2S (void) const
          /*
            \return number of transitions at 2-fold, transitional degenerate sites
          */
          {
            return p2S;
          }
          double P2V (void) const
          /*
            \return number of transitions at  2-fold, transversional degenerate sites
          */
          {
            return p2V;
          }
          double P4 (void) const
          /*
            \return number of transitions at 4-fold degenerate sites
          */
          {
            return p4;
          }
          double Q0 (void) const
          /*
            \return number of transversion at nondegenerate sites
          */
          {
            return q0;
          }
          double Q2S (void) const
          /*
            \return number of transversion at 2-fold, transitional degenerate sites
          */
          {
            return q2S;
          }
          double Q2V (void) const
          /*
            \return number of transversion at 2-fold, transversional sites
          */
          {
            return q2V;
          }
          double Q4 (void) const
          /*
            \return number of transversion at 4-fold degenerate sites
          */
          {
            return q4;
          }
          
          
          
            double L0 (void) const
          /*
            \return the number of nondegenerate sites compared
          */
          {
            return sites->L0();
          }
          double L2S (void) const
          /*
            \return the number of twofold, transitional-degenerate sites compared
          */
          {
            return sites->L2S();
          }
          double L2V (void) const
          /*
            \return the number of twofold, transversional-degenerate sites compared
          */
          {
            return sites->L2V();
          }
          double L4 (void) const
          /*
            \return the number of 4-fold degenerate sites compared
          */
          {
            return sites->L4();
          }
        
          
  double as (void) const
  /*
    \return corrected synonymous divergence at transitional-degenerate sites
  */
  {
    if (!isfinite (As))
      return 999.;
    return As;
  }
  double aa (void) const
  /*
    \return corrected nonsynonymous divergence at tranversioal- and non- degenerate sites
  */
  {
    if (!isfinite (Aa))
      return 999.;
    return Aa;
  }
  double bs (void) const
  /*
    \return corrected synonymous divergence at transversional- and fourfold-  degenerate sites
  */
  {
    if (!isfinite (Bs))
      return 999.;
    return Bs;
  }
  double ba (void) const
  /*
    \return corrected nonsynonymous divergence at transitional- and non- degenerate sites
  */
  {
    if (!isfinite (Ba))
      return 999.;
    return Ba;
  }
};

