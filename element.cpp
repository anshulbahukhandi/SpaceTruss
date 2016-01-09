/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CElement class implementation.

*********************************************/
#include "Element.h"

CElement::CElement ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_nSN = 0;
    m_nEN = 0;
    m_fArea = 0.0f;
    m_fE = 0.0f;
	m_ThermalCoeff=0.0f;
}

CElement::~CElement ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElement::GetData (int& nSN, int& nEN, float& fA,
                        double& fE,double& Tc) const
// ---------------------------------------------------------------------------
// Function: gets element-related values
// Input:    variables to hold start node#, end node #, x/s area, E
// Output:   start node#, end node #, x/s area, E
// ---------------------------------------------------------------------------
{
    nSN = m_nSN;
    nEN = m_nEN;
    fA = m_fArea;
    fE = m_fE;
	Tc= m_ThermalCoeff;
}

void CElement::SetData (const int nSN, const int nEN,
                        const float fA, const double fE, const double Tc)
// ---------------------------------------------------------------------------
// Function: sets element-related values
// Input:    variables holding values for start node#, end node #, 
//           x/s area, E
// Output:   modified value for start node#, end node #, x/s area, E
// ---------------------------------------------------------------------------
{
    m_nSN = nSN;
    m_nEN = nEN;
    m_fArea = fA;
    m_fE = fE;
	m_ThermalCoeff=Tc;
}
