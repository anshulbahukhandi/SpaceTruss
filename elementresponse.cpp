/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CElementResponse class implementation.

*********************************************/
#include "ElementResponse.h"

CElementResponse::CElementResponse ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fStrain = 0.0f;
    m_fStress = 0.0f;
    m_fForce = 0.0f;

}

CElementResponse::~CElementResponse ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CElementResponse::GetData (float& fE,
                                float& fS, float& fF) const
// ---------------------------------------------------------------------------
// Function: gets the element strain, stress and force
// Input:    variables to hold the values
// Output:   element strain, stress and force values
// ---------------------------------------------------------------------------
{
    fE = m_fStrain;
    fS = m_fStress;
    fF = m_fForce;

}

void CElementResponse::SetData (const float fE,
                                const float fS, const float fF)
// ---------------------------------------------------------------------------
// Function: sgets the element strain, stress and force
// Input:    variables holding the values
// Output:   modified element strain, stress and force values
// ---------------------------------------------------------------------------
{
    m_fStrain = fE;
    m_fStress = fS;
    m_fForce = fF;

}