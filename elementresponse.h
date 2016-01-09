/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CElementResponse class definition.

*********************************************/
#ifndef __RAJAN_ELEMENTRESPONSE_H__
#define __RAJAN_ELEMENTRESPONSE_H__

class CElementResponse
{
    public:
        // ctor and dtor
        CElementResponse ();
        ~CElementResponse ();

        // helper functions
        void Print () const;

        // accessor function
        void GetData (float& fStrain, float& fStress, 
                      float& fForce) const;

        // modifier functions
        void SetData (const float fStrain, const float fStress,
                      const float fForce);

    private:
        float m_fStrain; // strain
        float m_fStress; // stress
        float m_fForce;	 // force
};

#endif