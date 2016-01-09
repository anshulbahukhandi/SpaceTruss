/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CElement class definition.

*********************************************/
#ifndef __RAJAN_ELEMENT_H__
#define __RAJAN_ELEMENT_H__

class CElement
{
    public:
        // ctor and dtor
        CElement ();
        ~CElement ();

        // accessor function
        void GetData (int& nSN, int& nEN, float& fArea, 
                      double& fE,double& fTc) const;

        // modifier functions
        void SetData (const int nSN, const int nEN, const float fArea,
                      const double fE, const double fTc);

    private:
        int m_nSN;		// start node
        int m_nEN;		// end node
        float m_fArea;	// x/s area
        double m_fE;		// modulus of elasticity
        double m_ThermalCoeff;   // Coefficient of thermal expansion
};

#endif