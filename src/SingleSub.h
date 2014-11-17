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

#include "RedundancyCom95.h"

class SingleSub{
        
 private:
      double q0i, q2Si, q2Vi, q4i, q0j, q2Sj, q2Vj, q4j, p0i, p2Si, p2Vi,
      p4i, p0j, p2Sj, p2Vj, p4j;
      
   
 void Calculate (const RedundancyCom95 * sitesObj, const std::string & codon1,
                        const std::string & codon2)
  /*
    count up mutations between the codons
  */
  {
    string type_string;
    int pos = 0, k = 0, type=0;
    double d;
    for (k = 0; k <= 2; ++k)
      {
        if (codon1[k] != codon2[k])
	  {
	    pos = k;
	  }
      }
    type_string = TsTv (nucToInt(codon1[pos]), nucToInt(codon2[pos]));
    if(type_string == "Ts") type = 1;
    if(type_string == "Tv") type = 2;
    
//    if( type == "Unknown" )
//      {
//	ostringstream o;
//	o << "SingleSub.cc: mutation between " << codon1 << " and " << codon2
//	  << " at position " << pos
//	  << " is neither a transition nor a transversion";
//	throw SeqException(o.str().c_str());
//      }
    assert(type==1 || type==2);

    switch (pos)
      {
      case 0: //mutation at first position
	switch (type)
	  {
	  case 1: //transition
	    d = sitesObj->FirstNon(codon1);
	    if( d < 1. )
	      {
		d = sitesObj->First2S(codon1);
		if(d > 0.)
		  {
		    p2Si += 1.;
		  }
		else
		  {
		    p2Vi += 1.;
		  }
	      }
	    else
	      {
		p0i += 1.;
	      }
	    d = sitesObj->FirstNon(codon2);
	    if( d < 1. )
	      {
		d = sitesObj->First2S(codon2);
		if(d > 0.)
		  {
		    p2Sj += 1.;
		  }
		else
		  {
		    p2Vj += 1.;
		  }
	      }
	    else
	      {
		p0j += 1.;
	      }
	    break;
	  case 2: //transversion 
	    d = sitesObj->FirstNon(codon1);
	    if( d < 1. )
	      {
		d = sitesObj->First2V(codon1);
		if(d > 0.)
		  {
		    q2Vi += 1.;
		  }
		else
		  {
		    q2Si += 1.;
		  }
	      }
	    else
	      {
		q0i += 1.;
	      }
	    d = sitesObj->FirstNon(codon2);
	    if( d < 1. )
	      {
		d = sitesObj->First2V(codon2);
		if(d > 0.)
		  {
		    q2Vj += 1.;
		  }
		else
		  {
		    q2Sj += 1.;
		  }
	      }
	    else
	      {
		q0j += 1.;
	      }
	    break;
	  }
	break;
      case 1: //mutation at second position
	switch (type)
	  {
	  case 1: //transition
	    p0i += 1.0;
	    p0j += 1.0;
	    break;
	  case 2: //transversion 
	    q0i += 1.0;
	    q0j += 1.0;
	    break;
	  }
	break;
      case 2: //mutation at third position
	switch (type)
	  {
	  case 1: //transition
	    d = sitesObj->ThirdNon (codon1);
	    if ( d < 1.)
	      {
		d = sitesObj->ThirdFour (codon1);
		if ( d < 1. )
		  {
		    d = sitesObj->Third2S (codon1);
		    if ( d > 0. )
		      {
			p2Si += 1.;
		      }
		    else
		      {
			p2Vi += 1.;
		      }
		  } 
		else
		  {
		    p4i += 1.;
		  }
	      }
	    else
	      {
		p0i += 1.;
	      }

	    d = sitesObj->ThirdNon (codon2);
	    if ( d < 1.)
	      {
		d = sitesObj->ThirdFour (codon2);
		if ( d < 1. )
		  {
		    d = sitesObj->Third2S (codon2);
		    if ( d > 0. )
		      {
			p2Sj += 1.;
		      }
		    else
		      {
			p2Vj += 1.;
		      }
		  } 
		else
		  {
		    p4j += 1.;
		  }
	      }
	    else
	      {
		p0j += 1.;
	      }
	    break;
	  case 2: //transversion 
	    d = sitesObj->ThirdNon (codon1);
	    if ( d < 1.)
	      {
		d = sitesObj->ThirdFour (codon1);
		if ( d < 1. )
		  {
		    d = sitesObj->Third2V (codon1);
		    if ( d > 0. )
		      {
			q2Vi += 1.;
		      }
		    else
		      {
			q2Si += 1.;
		      }
		  } 
		else
		  {
		    q4i += 1.;
		  }
	      }
	    else
	      {
		q0i += 1.;
	      }

	    d = sitesObj->ThirdNon (codon2);
	    if ( d < 1.)
	      {
		d = sitesObj->ThirdFour (codon2);
		if ( d < 1. )
		  {
		    d = sitesObj->Third2V (codon2);
		    if ( d > 0. )
		      {
			q2Vj += 1.;
		      }
		    else
		      {
			q2Sj += 1.;
		      }
		  } 
		else
		  {
		    q4j += 1.;
		  }
	      }
	    else
	      {
		q0j += 1.;
	      }
	    break;
	  }
	break;
      }
      
  }
  
  public:      
//   SingleSub(const RedundancyCom95 * sitesObj,
//                	      const std::string & cod1,
//			      const std::string & cod2)
 void operator()(const RedundancyCom95 * sitesObj,
                	      const std::string & cod1,
			      const std::string & cod2)
  /*
    \param sitesObj an object of type Sequence::RedundancyCom95
    \param cod1 a std::string of length 3 representing a codon
    \param cod2 a std::string of length 3 representing a codon
    \note cod1 and cod2 lengths are verified by assert()
  */
  {
    assert(cod1.length() == 3 && cod2.length() == 3);
    q0i = q2Si = q2Vi = q4i = q0j = q2Sj = q2Vj = q4j = p0i = p2Si =
      p2Vi = p4i = p0j = p2Sj = p2Vj = p4j = 0.0;
    Calculate (sitesObj, cod1, cod2);
    
  }
   double P0(void) const
  /*
    \return number of transitions at non-degenerate sites in the codon
  */
  {
    return (p0i+p0j)/2.0;
  }


  double P2S(void) const
  /*
    \return number of transitions at transitional-degenerate sites in the codon
  */
  {
    return (p2Si+p2Sj)/2.0;
  }


  double P2V(void) const
  /*
    \return number of transitions at transversional-degenerate sites in the codon
  */
  {
    return (p2Vi+p2Vj)/2.0;
  }


  double P4(void) const
  /*
    \return number of transitions at fourfold-degenerate sites in the codon
  */
  {
    return (p4i+p4j)/2.0;
  }

  double Q0(void) const
  /*
    \return number of transversions at non-degenerate sites in the codon
  */
  {
    return (q0i+q0j)/2.0;
  }


  double Q2S(void) const
  /*
    \return number of transversions at transitional-degenerate sites in the codon
  */
  {
    return (q2Si+q2Sj)/2.0;
  }


  double Q2V(void) const
  /*
    \return number of transversions at transversional-degenerate sites in the codon
  */
  {
    return (q2Vi+q2Vj)/2.0;
  }


  double Q4(void) const
  /*
    \return number of transversions at fourfold-degenerate sites in the codon
  */
  {
    return (q4i+q4j)/2.0;
  }
};
