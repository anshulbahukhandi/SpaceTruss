/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CNode class definition.

*********************************************/
#ifndef __RAJAN_NODE_H__
#define __RAJAN_NODE_H__

class CNode
{
    public:
        // ctor and dtor
        CNode ();
        ~CNode ();

        // accessor functions
        void GetCoords (float&, float&, float&) const;    
        void GetFixity (int&, int&, int&) const;             
        void GetLoads (float&, float&, float&) const;    
		float GetTemperature()const;
		void GetNodalSettlement(float& , float& , float& )const;
        // modifier functions
        void SetCoords (const float, const float, const float);   
        void SetFixity (const int, const int, const int);    
        void SetLoads (const float, const float,const float);   
		void SetTemperature(const float);
		void SetNodalSettlement(const float , const float ,const float);
    private:
        float m_fXCoor;		// x-coordinate
        float m_fYCoor;		// y-coordinate
		float m_fZCoor;     // z-coordinate        
        int   m_nXFC;		// x-fixity code
        int   m_nYFC;		// y-fixity code
		int   m_nZFC;       // z-fixity code       
        float m_fXForce;	// x-force
        float m_fYForce;	// y-force
        float m_fZForce;    // z-force             
	    float m_fTemperature;   //Nodal Temperature
		float m_fXDisp;			// x nodal displacement/settlement
		float m_fYDisp;			// y nodal displacement/settlement
		float m_fZDisp;			//z nodal Displacement/settlement
};

#endif	