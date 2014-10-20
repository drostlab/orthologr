#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout


#include "Comeron95.h"

using namespace Rcpp;
using namespace std;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
// @export
// [[Rcpp::export]]
void gestimator(std::string file, std::string file_out = "", int maxHits = 3, 
                bool verbose = false, bool remove_all_gaps = false) {

   // read codon file
   ifstream in_stream (file.c_str());
   string line;
   vector<pair<string,string> > data;
   string sequence = "";
   pair<string,string> seq;
   while(getline(in_stream,line)){
           //in_stream >> line;
         
           size_t found = line.find(">");
           // new seq found
           if (found!=std::string::npos){
                // push pevious seq
                if(sequence != "") {
                        seq.second = sequence;
                        data.push_back(seq);
                        sequence = "";
                }
                // get name for next sequence
                seq.first = line.substr(1, line.find(" "));
           }
           else {
                // otherwise append to sequence
                sequence.append(line);
           }
   }
   seq.second = sequence;
   data.push_back(seq);


//  // check correct infile reading 
//   for(int i=0;i<data.size();i++){
//           Rcpp::Rcout << i << endl;
//           Rcpp::Rcout << data[i].first << endl;
//           Rcpp::Rcout << data[i].second << endl;
//           Rcpp::Rcout << data[i].second.length()<<endl;
//   }

   in_stream.close();
   

   // prepare output
   
   ofstream out_stream;
   if (file_out != "")
        out_stream.open(file_out.c_str());

   if (verbose == true)
    {
      if (file_out != "")
        {
          out_stream <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
          out_stream <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
        }
      else
        {
          Rcpp::Rcout <<"seq1\tseq2\tKa\tKs\tomega\tp0\tp2S\tp2V\tp4\tq0\tq2S\tq2V\tq4\t";
          Rcpp::Rcout <<"Aa\tAs\tBa\tBs\tL0\tL2S\tL2V\tL4"<<endl;
        }
     }
   
 
   // delete gaps
   if(remove_all_gaps){
           
   }
   
   // process ()
   
   RedundancyCom95 code_degeneracy;
   
   // foreach sequence pair in data
   for(int i=0;i<data.size();i++){
           for(int j=i+1;j<data.size();j++){
                   Comeron95 *C = new Comeron95(&data[i], &data[j], maxHits,&code_degeneracy);
                   if (file_out == "")
                        {
                          Rcpp::Rcout.precision(4); 
                          Rcpp::Rcout << data[i].first << '\t';
                          Rcpp::Rcout << data[j].first << '\t';
                          Rcpp::Rcout << C->dn() << '\t';
                          Rcpp::Rcout << C->ds() << '\t';
                          Rcpp::Rcout << C->ratio();
                          if (verbose)
                            {
                              Rcpp::Rcout <<'\t';
                              Rcpp::Rcout <<C->P0() << '\t' <<C->P2S()<<'\t';
                              Rcpp::Rcout <<C->P2V() <<'\t' <<C->P4()<<'\t';
                              Rcpp::Rcout <<C->Q0() << '\t' <<C->Q2S()<<'\t';
                              Rcpp::Rcout <<C->Q2V() << '\t' <<C->Q4()<<'\t';
                              Rcpp::Rcout <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
                              Rcpp::Rcout.precision(6);
                              Rcpp::Rcout <<C->L0() << '\t' <<C->L2S()<<'\t';
                              Rcpp::Rcout <<C->L2V() << '\t' <<C->L4()<<endl;
                            }
                          else
                            Rcpp::Rcout << endl;
                        }
                      else
                        {
                          out_stream.precision(4);
                          out_stream << data[i].first << '\t';
                          out_stream << data[j].first << '\t';
                          out_stream << C->dn() << '\t';
                          out_stream << C->ds() << '\t';
                          out_stream << C->ratio();
                          if (verbose)
                            {
                              out_stream<<'\t';
                              out_stream <<C->P0() << '\t' <<C->P2S()<<'\t';
                              out_stream <<C->P2V() <<'\t' <<C->P4()<<'\t';
                              out_stream <<C->Q0() << '\t' <<C->Q2S()<<'\t';
                              out_stream <<C->Q2V() << '\t' <<C->Q4()<<'\t';
                              out_stream <<C->aa()<<'\t'<<C->as()<<'\t'<<C->ba()<<'\t'<<C->bs()<<'\t';
                              out_stream.precision(6);
                              out_stream <<C->L0() << '\t' <<C->L2S()<<'\t';
                              out_stream <<C->L2V() << '\t' <<C->L4()<<endl;
                            }
                          else
                            out_stream<<endl;

                         }
                         
                        // manually deleting the previous results
                        //delete C;
                        C=NULL;
                 } // end j
                 
           } // end i
           
           out_stream.close();
} // end fun
  
   

   

   





/* IN: infile, outfile, verbose (if output should be written), maxhits, remove_all_gaps

analysis: gestimator
libsequence: comeron95, Alignment


gestimator.cc/Main:
- parse & set args -> args
- call gestimator.cc/process()
 
 gestimator.cc/process():
 In: pointer to filenames?, args, outstream
 - read infile
 - remove gaps with Alignment.cc/RemoveGaps()
 foreach pairs of input file:
        - call Comeron95.cc/Comeron95 Konstructor
        - read results from Comeron95-Object -> Names, ka, ks, ratio
        - print if verbose, print to file else
 
 Comeron95.cc/Comeron95 (Konstructor)
 In: 2 sequences, maxHits Value, ??/RedundancyCom95 pointer
 In optional: ??/GeneticsCode object, 2 weighting schemes
 - check input params (same length, max in 1-3)
 - assign default weighting scheme ??/GranthamWeights2/3(code), code without default
 - assign new sites object, init q0 to p4 with 0.0
 - call Comeron95/diverge(seqa,seqb,weights)
 - call ?/omega(seqa,seqb)
 
 Comeron95/diverge()
 In: 2 sequence Pointer, 2 weighting scheme pointer
 - loop overall codon pairs
        if(?/Different(codons))
                switch
                diff=1 : ?/SingleSub(sitesObj) -> Add predictedt changes to q0-p4
                diff=2 & maxhits >=2 : ?/TwoSubs(sitesObj, weight2) -> Add predictedt changes to q0-p4
                diff=3 & maxhits >2 : ?/ThreeSubs(sitesObj, weight3) -> Add predictedt changes to q0-p4
                if(some of q0-p4 is infinity) set 0.0
                
 Comeron95/omega()
 In: 2 seqs, from somewhere a site object
 - implements all equations given in the paper
 

*/