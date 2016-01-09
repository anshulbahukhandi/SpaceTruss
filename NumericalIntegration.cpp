/*********************************************
A  FEA  Template  Program
Copyright(c) 2004-06, S. D. Rajan
All rights reserved
*********************************************/

#include <cmath>
#include "NumericalIntegration.h"

//*****************************************************************
//************************** NATURAL COORDINATES ******************
//*****************************************************************
CGaussQuad::CGaussQuad ()
// ==================================================================
// Function: Constructor
//    Input: None
//   Output: None
// ==================================================================
{
    // location: order 1
    GC[0][0] =  0.0;

    // location: order 2
    GC[1][0] = -1.0/sqrt(3.0);
    GC[1][1] = -GC[1][0];

    // location: order 3
    GC[2][0] = -sqrt(3.0)/sqrt(5.0);
    GC[2][1] =  0.0;
    GC[2][2] = -GC[2][0];

    // location: order 4
    GC[3][0] = -sqrt(15.0+2.0*sqrt(30.0))/sqrt(35.0);
    GC[3][1] = -sqrt(15.0-2.0*sqrt(30.0))/sqrt(35.0);
    GC[3][2] = -GC[3][1];
    GC[3][3] = -GC[3][0];

    // location: order 5
    GC[4][0] = -sqrt(35.0+2.0*sqrt(70.0))/sqrt(63.0);
    GC[4][1] = -sqrt(35.0-2.0*sqrt(70.0))/sqrt(63.0);
    GC[4][2] = 0.0;
    GC[4][3] = -GC[4][1];
    GC[4][4] = -GC[4][0];

    // weight: order 1
    GW[0][0] = 2.0;

    // weight: order 2
    GW[1][0] = 1.0;
    GW[1][1] = 1.0;

    // weight: order 3
    GW[2][0] = 5.0/9.0;
    GW[2][1] = 8.0/9.0;
    GW[2][2] = GW[2][0];

    // weight: order 4
    GW[3][0] = 49.0/(6.0*(18.0+sqrt(30.0)));
    GW[3][1] = 49.0/(6.0*(18.0-sqrt(30.0)));
    GW[3][2] = GW[3][1];
    GW[3][3] = GW[3][0];

    // weight: order 5
    GW[4][0] = 5103.0/(50.0*(322.0+13.0*sqrt(70.0)));
    GW[4][1] = 5103.0/(50.0*(322.0-13.0*sqrt(70.0)));
    GW[4][2] = 128.0/225.0;
    GW[4][3] = GW[4][1];
    GW[4][4] = GW[4][0];
}

CGaussQuad::~CGaussQuad ()
// ==================================================================
// Function: Destructor
//    Input: None
//   Output: None
// ==================================================================
{
}

double CGaussQuad::GetLocation (int nOrder, int nPoint) const
// ==================================================================
// Function: Obtains the location for "nPoint" of "nOrder" rule
//    Input: nOrder and nPoint
//   Output: returns the location
// ==================================================================
{
    assert (nOrder >= 1 && nOrder <= MAXORDER);
    assert (nPoint >= 1 && nOrder <= MAXPOINT);
    return GC[nOrder-1][nPoint-1];
}

double CGaussQuad::GetWeight (int nOrder, int nPoint) const
// ==================================================================
// Function: Obtains the weight for "nPoint" of "nOrder" rule
//    Input: nOrder and nPoint
//   Output: returns the weight
// ==================================================================
{
    assert (nOrder >= 1 && nOrder <= MAXORDER);
    assert (nPoint >= 1 && nOrder <= MAXPOINT);
    return GW[nOrder-1][nPoint-1];
}

