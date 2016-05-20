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
#include <iostream>
#include <fstream>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout


#include "Comeron95.h"

using namespace Rcpp;
using namespace std;

// @export
// [[Rcpp::export]]
void gestimator(std::string file, std::string file_out = "", int maxHits = 3, 
                bool verbose = false, bool remove_all_gaps = false) {

   // read codon file -> GetData()
   ifstream in_stream (file.c_str());
   
   // check parameter
   if(!in_stream){
            Rcpp::stop("Input File not found!");
   }
   
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

		// if (Gapped(&data))
		bool gaps=false;
	    for(size_t i=0;i<data.size();i++){
			if( data[i].second.find("-") != std::string::npos){
				    gaps=true;
			}
	    }
		if(gaps){
		// RemoveGaps(data);

			/*
			  Modifies the data vector to remove all positions that
			  contain the gap character'-'.
			  \param data vector<T> to modify
			*/

		   size_t i, j,datasize=data.size();
		   int length = data[0].second.length ();
		   vector<string> ungapped_sequences(data.size());
		   bool site_is_gapped;

		  for (i = 0; i < length; ++i)
		    {        //iterate over sites
		      for ( j = 0, site_is_gapped = 0;
		            j < datasize;  ++j)
		        {
		          if (data[j].second[i] == '-')
		            {
		              site_is_gapped = 1;
		              j = datasize;
		            }
		        }
		      if (!(site_is_gapped))
		        {
		          for ( j = 0 ; j != data.size();  ++j)
		            ungapped_sequences[j] += data[j].second[i];
		        }
		    }

		  //redo the data
		  for (j = 0; j != datasize ; ++j)
		    {
		      data[j].second = ungapped_sequences[j];
		    }
        
		}
   }
   
   // process ()
   
   RedundancyCom95 code_degeneracy;
   
   // foreach sequence pair in data
   for(size_t i=0;i<data.size();i++){
           for(size_t j=i+1;j<data.size();j++){

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
