/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CTruss class definition.

*********************************************/
#ifndef __RAJAN_TRUSS_H__
#define __RAJAN_TRUSS_H__

#include <fstream>
#include <iostream>
#include "F:\ASU Courses\cee532\Library\Library\vectortemplate.h"
#include "F:\ASU Courses\cee532\Library\Library\matrixtemplate.h"
#include "node.h"
#include "element.h"
#include "nodalresponse.h"
#include "elementresponse.h"

class CTruss
{
    public:
        // ctor and dtor
        CTruss ();
        ~CTruss ();
        enum ERRORCODE {INVALIDNUMNODES, INVALIDNUMELEMENTS, INVALIDDEBUGCODE,
                        INVALIDNODENUM,  INVALIDNODALFIXITY, INVALIDELEMENTNUM,
                        INVALIDCSAREA,   INVALIDYM,          MISSINGEND,
                        UNSTABLETRUSS,   INVALIDINPUT,       INVALIDCOMMANDLINE, INDETERMINANT};		//CHANGED THIS

        void Banner (std::ostream& OF) const;
        void PrepareIO (int argc, char* argv[]);
        void ReadProblemSize ();
        void ReadTrussModel ();
        void Analyze ();
        void TerminateProgram ();
		void GetErrors( );
    
private:
        int m_nNodes;		// number of nodes
        int m_nElements;	// number of elements
        int m_nDOF;			// total degrees-of-freedom					 
        int m_nDebugLevel;	// debugging level
        int m_nLineNumber;	// current line number in input file         

        // data storage for 
        CVector<CNode> m_NodalData;                     // nodal data
        CVector<CElement> m_ElementData;                // element data
        CVector<CNodalResponse> m_NodalResponseData;    // nodal response
        CVector<CElementResponse> m_ElementResponseData;// element response

        std::ifstream m_FileInput;	// File Input
        std::ofstream m_FileOutput;	// File Output

        CMatrix<double> m_dSSM;	// structural stiffness matrix
        CMatrix<double> m_dSND;	// structural nodal displacements
        CMatrix<double> m_dSNF;	// structural nodal forces
        float m_fAE;			 //Absolute Error					CHANGED THIS
		float m_fRE;			//Relative Error					CHANGED THIS
        void SetSize ();
        void ConstructK ();
        void ImposeBC ();
        void Solve ();
        void Response ();
        void CreateOutput ();
        void SuppressDOF (const int);
        void ErrorHandler (ERRORCODE) const; // handles input-related errors
};

#endif