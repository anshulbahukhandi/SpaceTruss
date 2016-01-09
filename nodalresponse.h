/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CNodalResponse class definition.

*********************************************/
#ifndef __RAJAN_NODALRESPONSE_H__
#define __RAJAN_NODALRESPONSE_H__

class CNodalResponse
{
    public:
        // ctor and dtor
        CNodalResponse ();
        ~CNodalResponse ();

        // helper functions
        void Print () const;

        // accessor function
        void GetDisplacements (float& fXDisp, float& fYDisp,float& fZDisp) const;

        // modifier functions
        void SetDisplacements (const float fXDisp, const float fYDisp,const float fZDisp);
		void CNodalResponse::GetSupportReaction(double& fXSR,double& fYSR,double& fZSR);						//CHANGED THIS 
		void CNodalResponse::SetSupportReaction(const double fXSR,const double fYSR,const double fZSR);			//CHANGED THIS

    private:
        float m_fXDisp;	// x-displacement
        float m_fYDisp;	// y-displacement
		float m_fZDisp; //z -displacement
		float m_fXSupportReaction;   //support reaction in X Direction					//CHANGED THIS
		float m_fYSupportReaction;   //support reaction in Y Direction					//CHANGED THIS
		float m_fZSupportReaction;	//support reaction in Z Direction					//CHANGED THIS
};

#endif