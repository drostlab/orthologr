#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

#include "TwoSubs.h"

class ThreeSubs{
        private:
      double p0, p2S, p2V, p4, q0, q2S, q2V, q4;

  void
  Calculate (const RedundancyCom95 * sitesObj,
                        const std::string *intermediates,
                        const std::string & codon1, const std::string & codon2,
                        double w_path1, double w_path2, double w_path3,
                        double w_path4, double w_path5, double w_path6)
  /*
    calculate number of changes along each pathway
  */
  {
    //arrays to store transition / transversion values per branch
    double p0_b[15], p2S_b[15], p2V_b[15], p4_b[15];
    double q0_b[15], q2S_b[15], q2V_b[15], q4_b[15];

    //initialuze the SingleSub class for the first path
    SingleSub Single;
    //there are 15 branches to iterate over--you'll need to draw a picture...
    //there is a picture in the documentation
    for (int i = 0; i <= 14; ++i)
      {
        switch (i)
          {
          case 0:
            //already initialized above
            // *Single =  SingleSub (sitesObj, codon1,intermediates[0]);
            Single(sitesObj, codon1,intermediates[0]);
            break;
          case 1:
            Single (sitesObj, intermediates[0],
                    intermediates[1]);
            break;
          case 2:
            Single (sitesObj, intermediates[1],
                    codon2);
            break;
          case 3:
            Single (sitesObj, intermediates[0],
                    intermediates[2]);
            break;
          case 4:
            Single (sitesObj, intermediates[2],
                    codon2);
            break;
          case 5:
            Single (sitesObj, codon1,
                    intermediates[3]);
            break;
          case 6:
            Single (sitesObj, intermediates[3],
                    intermediates[4]);
            break;
          case 7:
            Single (sitesObj, intermediates[4],
                    codon2);
            break;
          case 8:
            Single (sitesObj, intermediates[3],
                    intermediates[5]);
            break;
          case 9:
            Single (sitesObj, intermediates[5],
                    codon2);
            break;
          case 10:
            Single (sitesObj, codon1,
                    intermediates[6]);
            break;
          case 11:
            Single (sitesObj, intermediates[6],
                    intermediates[7]);
            break;
          case 12:
            Single (sitesObj, intermediates[7],
                    codon2);
            break;
          case 13:
            Single (sitesObj, intermediates[6],
                    intermediates[8]);
            break;
          case 14:
            Single (sitesObj, intermediates[8],
                    codon2);
            break;
          }
        p0_b[i] = Single.P0();
        p2S_b[i] =Single.P2S();
        p2V_b[i] =Single.P2V();
        p4_b[i] = Single.P4();
        q0_b[i] = Single.Q0();
        q2S_b[i] =Single.Q2S();
        q2V_b[i] =Single.Q2V();
        q4_b[i] = Single.Q4();
      }
    //sum up changes along each branch, weighting by the
    //weight factor for each path
    p0 = (p0_b[0] + p0_b[1] + p0_b[2]) * w_path1
         + (p0_b[0] + p0_b[3] +  p0_b[4]) * w_path2
         + (p0_b[5] + p0_b[6] + p0_b[7]) * w_path3
         + (p0_b[5] + p0_b[8] + p0_b[9]) * w_path4
         + (p0_b[10] + p0_b[11] + p0_b[12]) * w_path5
         + (p0_b[10] + p0_b[13] + p0_b[14]) * w_path6;

    p2S = (p2S_b[0] + p2S_b[1] + p2S_b[2]) * w_path1
          + (p2S_b[0] + p2S_b[3] + p2S_b[4]) * w_path2
          + (p2S_b[5] + p2S_b[6] + p2S_b[7]) * w_path3
          + (p2S_b[5] + p2S_b[8] + p2S_b[9]) * w_path4
          + (p2S_b[10] + p2S_b[11] + p2S_b[12]) * w_path5
          + (p2S_b[10] + p2S_b[13] + p2S_b[14]) * w_path6;

    p2V = (p2V_b[0] + p2V_b[1] + p2V_b[2]) * w_path1
          + (p2V_b[0] + p2V_b[3] + p2V_b[4]) *   w_path2
          + (p2V_b[5] + p2V_b[6] + p2V_b[7]) * w_path3
          + (p2V_b[5] + p2V_b[8] + p2V_b[9]) * w_path4
          + (p2V_b[10] + p2V_b[11] + p2V_b[12]) * w_path5
          + (p2V_b[10] + p2V_b[13] + p2V_b[14]) * w_path6;

    p4 = (p4_b[0] + p4_b[1] + p4_b[2]) * w_path1
         + (p4_b[0] + p4_b[3] +  p4_b[4]) * w_path2
         + (p4_b[5] + p4_b[6] + p4_b[7]) * w_path3
         + (p4_b[5] + p4_b[8] +  p4_b[9]) * w_path4
         + (p4_b[10] + p4_b[11] + p4_b[12]) * w_path5
         + (p4_b[10] + p4_b[13] + p4_b[14]) * w_path6;

    q0 = (q0_b[0] + q0_b[1] + q0_b[2]) * w_path1
         + (q0_b[0] + q0_b[3] +  q0_b[4]) * w_path2
         + (q0_b[5] + q0_b[6] + q0_b[7]) * w_path3
         + (q0_b[5] + q0_b[8] + q0_b[9]) *  w_path4
         + (q0_b[10] + q0_b[11] + q0_b[12]) * w_path5
         +  (q0_b[10] + q0_b[13] + q0_b[14]) * w_path6;

    q2S = (q2S_b[0] + q2S_b[1] + q2S_b[2]) * w_path1
          + (q2S_b[0] + q2S_b[3] + q2S_b[4]) *   w_path2
          + (q2S_b[5] + q2S_b[6] + q2S_b[7]) * w_path3
          + (q2S_b[5] + q2S_b[8] + q2S_b[9]) * w_path4
          + (q2S_b[10] +q2S_b[11] + q2S_b[12]) * w_path5
          + (q2S_b[10] + q2S_b[13] + q2S_b[14]) * w_path6;

    q2V = (q2V_b[0] + q2V_b[1] + q2V_b[2]) * w_path1
          + (q2V_b[0] + q2V_b[3] +        q2V_b[4]) *  w_path2
          + (q2V_b[5] + q2V_b[6] + q2V_b[7]) * w_path3
          + (q2V_b[5] + q2V_b[8] + q2V_b[9]) * w_path4
          + (q2V_b[10] +  q2V_b[11] + q2V_b[12]) * w_path5
          + (q2V_b[10] + q2V_b[13] + q2V_b[14]) * w_path6;

    q4 = (q4_b[0] + q4_b[1] + q4_b[2]) * w_path1
         + (q4_b[0] + q4_b[3] +  q4_b[4]) * w_path2
         + (q4_b[5] + q4_b[6] + q4_b[7]) * w_path3
         + (q4_b[5] + q4_b[8] + q4_b[9]) *  w_path4
         + (q4_b[10] + q4_b[11] + q4_b[12]) * w_path5
         +  (q4_b[10] + q4_b[13] + q4_b[14]) * w_path6;
  }
  
        public:
      void operator()(const RedundancyCom95 * sitesObj,
                              const std::string &codon1, const std::string &codon2,
                               GranthamWeights3 *weights3)
  /*
    \param sitesObj an object of type Sequence::RedundancyCom95
    \param code see Sequence::GeneticCodes for valid values
    \param codon1 a std::string of length 3
    \param codon2 a std::string of length 3
    \param weights3 a weighting scheme for the pathways
    \note length of codons is checked by assert()
   */
  {
      
               
    assert(codon1.length() == 3 && codon2.length() == 3);
    string intermediates[9];
    weights3->Intermediates3(intermediates,codon1,codon2);
    p0 = p2S = p2V = p4 = q0 = q2S = q2V = q4 = 0.0;

    double weights[6];
    weights3->Calculate(weights, codon1,codon2);
   // double *weights = weights3->weights();
    Calculate (sitesObj, intermediates, codon1, codon2, weights[0],
               weights[1], weights[2], weights[3], weights[4], weights[5]);
               
//       Rcpp::Rcout << "q0  "<<q0<<endl;
//      Rcpp::Rcout << "q2S  "<<q2S<<endl;
//      Rcpp::Rcout << "q2V  " << q2V<<endl;
//      Rcpp::Rcout << "q4  " << q4 <<endl;
//  
//      Rcpp::Rcout << "p0  " << p0<<endl;
//      Rcpp::Rcout << "p2S  " << p2S<<endl;
//      Rcpp::Rcout << "p2V  " << p2V<<endl;
//      Rcpp::Rcout << "p4  " << p4<<endl;
  }


        double
      P0 (void) const
      /*
        \return number of transitions at non-degenerate sites in the codon
      */
      {
        return p0;
      }

      double
      P2S (void) const
      /*
        \return number of transitions at transitional-degenerate sites in the codon
      */
      {
        return p2S;
      }

      double
      P2V (void) const
      /*
        \return number of transitions at transversional-degenerate sites in the codon
      */
      {
        return p2V;
      }

      double
      P4 (void) const
      /*
        \return number of transitions at fourfold-degenerate sites in the codon
      */
      {
        return p4;
      }

      double
      Q0 (void) const
      /*
        \return number of transversions at non-degenerate sites in the codon
      */
      {
        return q0;
      }

      double
      Q2S (void) const
      /*
        \return number of transversions at transitional-degenerate sites in the codon
      */
      {
        return q2S;
      }

      double
      Q2V (void) const
      /*
        \return number of transversions at transversional-degenerate sites in the codon
      */
      {
        return q2V;
      }

      double
      Q4 (void) const
      /*
        \return number of transversions at fourfold-degenerate sites in the codon
      */
      {
        return q4;
      }
};
