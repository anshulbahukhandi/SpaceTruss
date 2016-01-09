/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-06, S. D. Rajan
All rights reserved
*********************************************/
#ifndef __RAJAN_GAUSSQUAD__H__
#define __RAJAN_GAUSSQUAD__H__

#include <cassert>

const int MAXORDER = 5;
const int MAXPOINT = 5;

class CGaussQuad
{
	public:
        CGaussQuad ();
        ~CGaussQuad ();

        // accessor functions
        double GetWeight   (int nOrder, int nPoint) const;
        double GetLocation (int nOrder, int nPoint) const;

    private:
        double GC[MAXORDER][MAXPOINT];
        double GW[MAXORDER][MAXPOINT];
};

#endif
