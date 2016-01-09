/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CNode class implementation.

*********************************************/
#include "Node.h"

CNode::CNode ()
// ---------------------------------------------------------------------------
// Function: default constructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
    m_fXCoor=0.0f;
    m_fYCoor=0.0f;
	m_fZCoor=0.0f;
    m_nXFC=0;
    m_nYFC=0;
	m_nZFC=0;
    m_fXForce=0.0f;
    m_fYForce=0.0f;
	m_fZForce=0.0f;
	m_fTemperature=0.0f;
	m_fXDisp=0.0f;
	m_fYDisp=0.0f;
	m_fZDisp=0.0f;
}

CNode::~CNode ()
// ---------------------------------------------------------------------------
// Function: destructor
// Input:    none
// Output:   none
// ---------------------------------------------------------------------------
{
}

void CNode::GetCoords (float& fX, float& fY,float& fZ) const
// ---------------------------------------------------------------------------
// Function: gets nodal coordinates
// Input:    variables to hold x-coordinate, y-coordinate values
// Output:   x-coordinate, y-coordinate values
// ---------------------------------------------------------------------------
{
    fX = m_fXCoor;
    fY = m_fYCoor;
	fZ=m_fZCoor;
}
float CNode::GetTemperature()const
{
	return m_fTemperature;
}

void CNode::SetTemperature(const float temp)
{
	m_fTemperature=temp;
}

void CNode::GetFixity (int& nXFC, int& nYFC, int& nZFC) const
// ---------------------------------------------------------------------------
// Function: gets nodal fixity conditions
// Input:    variables to hold x-fixity, y-fixity values
// Output:   x-fixity, y-fixity values
// ---------------------------------------------------------------------------
{
    nXFC = m_nXFC;
    nYFC = m_nYFC;
	nZFC = m_nZFC;
}

void CNode::GetLoads (float& fXF, float& fYF,float& fZF) const
// ---------------------------------------------------------------------------
// Function: gets nodal loads
// Input:    variables to hold x-force, y-force values
// Output:   x-force, y-force values
// ---------------------------------------------------------------------------
{
    fXF = m_fXForce;
    fYF = m_fYForce;
	fZF = m_fZForce;
}
void CNode::GetNodalSettlement(float& fXD, float& fYD, float& fZD)const
{
	fXD=m_fXDisp;
	fYD=m_fYDisp;
	fZD=m_fZDisp;
}

void CNode::SetNodalSettlement(const float fXD, const float fYD,const float fZD)
{
	m_fXDisp=fXD;
	m_fYDisp=fYD;
	m_fZDisp=fZD;
}
void CNode::SetCoords (const float fX, const float fY,const float fZ)
// ---------------------------------------------------------------------------
// Function: sets nodal coordinates
// Input:    variables holding x-coordinate, y-coordinate values
// Output:   modified x-coordinate, y-coordinate values
// ---------------------------------------------------------------------------
{
    m_fXCoor = fX;
    m_fYCoor = fY;
	m_fZCoor = fZ;
}

void CNode::SetFixity (const int nXFC, const int nYFC,const int nZFC)
// ---------------------------------------------------------------------------
// Function: sets nodal fixities
// Input:    variables x-fixity, y-fixity values
// Output:   modified x-fixity, y-fixity values
// ---------------------------------------------------------------------------
{
    m_nXFC = nXFC;
    m_nYFC = nYFC;
	m_nZFC = nZFC;
}

void CNode::SetLoads (const float fXF, const float fYF,const float fZF)
// ---------------------------------------------------------------------------
// Function: sets nodal loads
// Input:    variables holding x-force, y-force values
// Output:   modified x-force, y-force values
// ---------------------------------------------------------------------------
{
    m_fXForce = fXF;
    m_fYForce = fYF;
	m_fZForce = fZF;
}

