/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis

Contains CTruss class implementation.

*********************************************/
#include <iomanip>
#include<iostream>
#include <sstream>
#include "truss.h"
#include "F:\ASU Courses\cee532\Space Truss\Space Truss\MatToolBox.h"
#include "F:\ASU Courses\cee532\Library\library\fileio.h"
#include "F:\ASU Courses\cee532\Library\Library\parser.h"
const int MAXCHARS = 80;
std::string szInputString;
std::string szComment ("**");                          
CMatToolBox<double> MTB;

/* ==================================================================
======================= CTruss class =============================
================================================================== */

CTruss::CTruss ()
	// ---------------------------------------------------------------------------
	// Function: default ctor
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	m_nNodes = 0;
	m_nElements = 0;
	m_nDOF = 0;
	m_nDebugLevel = 0;
	m_nLineNumber = 0;                                           
}

CTruss::~CTruss ()
	// ---------------------------------------------------------------------------
	// Function: destructor
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
}

void CTruss::Banner (std::ostream& OF) const
	// ---------------------------------------------------------------------------
	// Function: prints the program banner on the output stream                     
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	OF << '\n';
	OF << "\t\t--------------------------------------------" << '\n';
	OF << "\t\t         Planar Truss Analysis Program      " << '\n';
	OF << "\t\tIntroduction to Structural Analysis & Design" << '\n';
	OF << "\t\t           (c) 2000-13, S. D. Rajan         " << '\n';
	OF << "\t\t         Enhanced By: Anshul Bahukhandi        " << '\n';
	OF << "\t\t--------------------------------------------" << '\n';
}

void CTruss::PrepareIO (int argc, char* argv[])
	// ---------------------------------------------------------------------------
	// Function: opens the input and output files
	// Input:    command line arguments (currently unused)
	// Output:   none
	// ---------------------------------------------------------------------------
{
	if (argc == 1)
	{   
		// open the input file
		OpenInputFileByName ("Complete input file name: ", m_FileInput,           
			std::ios::in);

		// open the output file
		OpenOutputFileByName ("Complete output file name: ", m_FileOutput,
			std::ios::out);


	}
	else
	{
		ErrorHandler (INVALIDCOMMANDLINE);
	}

	// print banner
	Banner (m_FileOutput);
}

void CTruss::ReadProblemSize ()                                 
	// ---------------------------------------------------------------------------
	// Function: Reads the size of the problem being solved
	// Input:    none

	// Output:   none
	// ---------------------------------------------------------------------------
{
	CParser Parse;

	// header line
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,      
		MAXCHARS, szComment))
		ErrorHandler (INVALIDINPUT);
	// read the problem description
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,     
		MAXCHARS, szComment))
		ErrorHandler (INVALIDINPUT);

	// header line
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,    
		MAXCHARS, szComment))
		ErrorHandler (INVALIDINPUT);						
	// read number of nodes, elements and debug level
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,     
		MAXCHARS, szComment))    
		ErrorHandler (INVALIDINPUT);
	std::istringstream szFormatString (szInputString);
	szFormatString >> m_nNodes >> m_nElements >> m_nDebugLevel;
	if (szFormatString.fail() || szFormatString.bad())
		ErrorHandler (INVALIDINPUT);

	// check data for validity
	if (m_nNodes <= 1) ErrorHandler (INVALIDNUMNODES);
	if (m_nElements <= 0) ErrorHandler (INVALIDNUMELEMENTS);
	if (m_nDebugLevel < 0 || m_nDebugLevel > 1) ErrorHandler (INVALIDDEBUGCODE);

	// dynamic memory allocations for arrays
	SetSize ();
}

void CTruss::SetSize ()
	// ---------------------------------------------------------------------------
	// Function: Carries out memory allocation for all the major arrays
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	// allocate space for nodal data
	m_NodalData.SetSize (m_nNodes);
	// allocate space for nodal response data
	m_NodalResponseData.SetSize (m_nNodes);
	// allocate space for element data
	m_ElementData.SetSize (m_nElements);
	// allocate space for element response data
	m_ElementResponseData.SetSize (m_nElements);

	// allocate and initialize major matrices
	m_nDOF = 3*m_nNodes;                       
	m_dSSM.SetSize (m_nDOF, m_nDOF);
	m_dSND.SetSize (m_nDOF, 1);
	m_dSNF.SetSize (m_nDOF, 1);
	m_dSSM.Set (0.0);
	m_dSND.Set (0.0);
	m_dSNF.Set (0.0);
}

void CTruss::ReadTrussModel ()
	// ---------------------------------------------------------------------------
	// Function: Read the truss model data from the input file
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	// read problem size
	ReadProblemSize ();

	int i, nN;
	float fX, fY,fZ;       
	CParser Parse;

	// header line                                       
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,       
		MAXCHARS, szComment))
		ErrorHandler (INVALIDINPUT);
	// read nodal coordinates
	for (i=1; i <= m_nNodes; i++)
	{
		if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString,   
			MAXCHARS, szComment))
			ErrorHandler (INVALIDINPUT);
		else
		{
			std::istringstream szFormatString (szInputString);                 
			szFormatString >> nN >> fX >> fY>> fZ;				
			if (szFormatString.fail() || szFormatString.bad())
				ErrorHandler (INVALIDINPUT);
		}
		if (nN <= 0 || nN > m_nNodes) 
			ErrorHandler (INVALIDNODENUM);
		m_NodalData(nN).SetCoords (fX, fY,fZ);

	}

	// header line                  
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
		MAXCHARS, szComment))
		ErrorHandler (INVALIDINPUT);
	// read nodal fixity conditions and Nodal Displacements 
	float fXD=0.0f;
	float fYD=0.0f;
	float fZD=0.0f;     
	int nXFC=0;
	int nYFC=0;
	int nZFC=0;    
	for (i=1; i <= m_nNodes+1; i++)
	{
		if (Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
			MAXCHARS, szComment))
		{
 			if(szInputString=="*nodal loads")												  //CHANGED THIS
				break;

			else
			{ 
				std::istringstream szFormatString (szInputString);
				szFormatString >> nN >> nXFC >> nYFC >>nZFC>> fXD >> fYD>> fZD;    
				if (szFormatString.fail() || szFormatString.bad())
					ErrorHandler (INVALIDINPUT);
				m_NodalData(nN).SetFixity (nXFC, nYFC, nZFC);   
				m_NodalData(nN).SetNodalSettlement(fXD , fYD , fZD);
			}

		}
		/*int nX,nY,nZ,TF=0;											was trying this condition but failed
		for(int i=1;i<=m_nNodes;i++)
			{
				m_NodalData(i).GetFixity(nX,nY,nZ);
				if(nX==1) TF++;
				if(nY==1) TF++;
				if(nZ==1) TF++;
		}
		if((m_nElements+TF)>3*m_nNodes)
		{
			ErrorHandler(INDETERMINANT);
		}*/

		if (nN <= 0 || nN > m_nNodes) ErrorHandler (INVALIDNODENUM);
		if (nXFC < 0 || nXFC > 1) ErrorHandler (INVALIDNODALFIXITY);
		if (nYFC < 0 || nYFC > 1) ErrorHandler (INVALIDNODALFIXITY);
		if (nZFC < 0 || nZFC > 1) ErrorHandler (INVALIDNODALFIXITY);

	}
	// read nodal forces and Nodal temperature if any
	float fXF=0.0;
	float fYF=0.0f;
	float fZF=0.0f;
	float fTemp=0.0f;         
	for (i=1; i <= m_nNodes+1; i++)
	{

		if (Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
			MAXCHARS, szComment))
		{
			if(szInputString=="*element data") 	break;            //CHANGED THIS

			else
			{ std::istringstream szFormatString (szInputString);
			szFormatString >> nN >> fXF >> fYF>> fZF>>fTemp;
			m_NodalData(nN).SetTemperature(fTemp);
			m_NodalData(nN).SetLoads(fXF,fYF,fZF);
			}
		}
	}

																	//CHANGED THIS ..REMOVED ONE READNEXTLINE function			
	// read element data  
	int nE, nSN, nEN;
	float fA;
	double ThermalCoefficient,dE;
	for (i=1; i <= m_nElements; i++)
	{
		if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
			MAXCHARS, szComment))
			ErrorHandler (INVALIDINPUT);
		else
		{
			std::istringstream szFormatString (szInputString);
			szFormatString >> nE >> nSN >> nEN >> fA >> dE>>ThermalCoefficient;
			if (szFormatString.fail() || szFormatString.bad())
				ErrorHandler (INVALIDINPUT);
		}
		if (nE <= 0 || nE > m_nElements) ErrorHandler (INVALIDELEMENTNUM);
		if (nSN <= 0 || nSN > m_nNodes) ErrorHandler (INVALIDNODENUM);
		if (nEN <= 0 || nEN > m_nNodes) ErrorHandler (INVALIDNODENUM);
		if (fA <= 0.0f) ErrorHandler (INVALIDCSAREA);
		if (dE <= 0.0f) ErrorHandler (INVALIDYM);
		m_ElementData(nE).SetData (nSN, nEN, fA, dE,ThermalCoefficient);
	}

	// end keyword?
	if (!Parse.ReadNextLine (m_FileInput, m_nLineNumber, szInputString, 
		MAXCHARS, szComment))
		ErrorHandler ((INVALIDINPUT));
	if (szInputString.substr(0,4) != "*end")
		ErrorHandler (MISSINGEND);

	// construct structural nodal load vector
	for (i=1; i <= m_nNodes; i++)
	{
		m_NodalData(i).GetLoads (fXF, fYF,fZF);
		m_dSNF(3*i-2,1) = fXF; 
		m_dSNF(3*i-1,1) = fYF;
		m_dSNF(3*i,1) = fZF;
	}
}

void CTruss::Analyze ()
	// ---------------------------------------------------------------------------
	// Function: Implements the FEA steps
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	// construct structural stiffness matrix
	ConstructK ();

	// impose boundary conditions
	ImposeBC ();

	// solve for the nodal displacements
	Solve ();

	// compute element response
	Response ();

	// create output file
	CreateOutput ();
}

void CTruss::ConstructK ()                              
	// ---------------------------------------------------------------------------
	// Function: Constructs the structural stiffness matrix and thermal load matrix 
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	int i, j, k;
	CVector<int> nE(6);
	CMatrix<double> dkl(2,2), dkg(6,6);
	CMatrix<double> dT(2,6), dTT(6,2);
	CMatrix<double> dTemp(2,6);
	float fX1, fX2, fY1, fY2,fZ1,fZ2, fL;
	float fA,TempDiff;
	int nSN, nEN;
	double da,fE;
	// initialize
	dT.Set (0.0);
	dTT.Set (0.0);
	CMatrix<double>dSNTF(m_nDOF,1);
	dSNTF.Set(0.0);	

	// loop thro' all elements
	for (i=1; i <= m_nElements; i++)
	{   
		m_ElementData(i).GetData(nSN,nEN,fA,fE,da);
		m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
		m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
		fL = sqrt((fX2-fX1)*(fX2-fX1) +
			(fY2-fY1)*(fY2-fY1)+(fZ2-fZ1)*(fZ2-fZ1));
		// form AE/L 
		double dAEOL = double(fA*fE/fL);       
		// construct klocal
		dkl(1,1) = dAEOL;
		dkl(1,2) = -dAEOL;
		dkl(2,1) = dkl(1,2);
		dkl(2,2) = dkl(1,1);
		// local-to-global transformation matrix
		double dl = double((fX2-fX1)/fL);
		double dm = double((fY2-fY1)/fL);
		double dn = double((fZ2-fZ1)/fL);
		dT(1,1) = dl; dT(1,2) = dm; dT(1,3)=dn;
		dT(2,4) = dl; dT(2,5) = dm; dT(2,6)=dn;

		//....................Constructing the global thermal force matrix......................
		nE(1) = 3*nSN-2; nE(2) = 3*nSN-1; nE(3) = 3*nSN;
		nE(4) = 3*nEN-2; nE(5) = 3*nEN-1; nE(6) = 3*nEN;


		TempDiff=m_NodalData(nEN).GetTemperature()-m_NodalData(nSN).GetTemperature();
		CMatrix<double>dEF (2,1);   //Thermal force in elements
		CMatrix<double>dNTF (6,1);
		dEF(1,1)=-fA*fE*da*TempDiff;
		dEF(2,1)=fA*fE*da*TempDiff;

		CMatrix<double>dSNTF(m_nDOF,1);
		dSNTF.Set(0.0f);

		//Adding the Thermal forces to the given structural nodal force matrix
		if(MTB.Transpose(dT,dTT))
		{

			if(MTB.Multiply(dTT,dEF,dNTF))      //constructing global thermal matrix
			{   
				for (int l=1;l<=6;l++)
				{ int nrow=nE(l);

				dSNTF(nrow,1)=dNTF(l,1);
				}


			}
		}
		for (int l=1; l<=m_nDOF ; l++)    //adding  the thermal matrix to nodal force matrix
		{
			m_dSNF(l,1)+=dSNTF(l,1);
		}
// construct k(6x6) in two steps   .........................................


		// transpose of the T matrix
		if(MTB.Transpose (dT, dTT))
		{                     
			// form k'*T
			if(MTB.Multiply ( dkl, dT,dTemp))                    
				// construct kglobal=T(T)*k'*T              
			{ 
				if(MTB.Multiply ( dTT, dTemp,dkg))
				{                                            

					// assemble into structural K

					for (j=1; j<= 6; j++)
					{
						int nRow = nE(j);
						for (k=1; k <= 6; k++)
						{
							int nCol = nE(k);
							m_dSSM(nRow, nCol) += dkg(j,k);

						}
					}
				}
			}
		}

		// debug?
		if (m_nDebugLevel == 1)
		{
			std::ostringstream szPrompt;
			szPrompt << "Stiffness Matrix for element " << i;
			MTB. PrintMatrixRowWise (dkg, szPrompt.str(), m_FileOutput);  //why szprompt.str() here and not below???!!
		}
	}
// debug?
	if (m_nDebugLevel == 1)
	{
		MTB.PrintMatrixRowWise (m_dSSM, "Structural Stiffness (Before BCs)",
			m_FileOutput);                    
		std::cout<<"4";
	}
}

void CTruss::ImposeBC ()
	// ---------------------------------------------------------------------------
	// Function: Imposes the essential boundary conditions
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	int i;
	int nXFC, nYFC,nZFC;
	float fXD=0.0f;
	float fYD=0.0f;
	float fZD=0.0f;
	// loop thro' all nodes
	for (i=1; i <= m_nNodes; i++)
	{

		//Imposing boundary conditions on nodal displacements


		m_NodalData(i).GetNodalSettlement(fXD,fYD,fZD); 

		m_NodalData(i).GetFixity (nXFC, nYFC,nZFC);
		if (nXFC == 1)
		{
			int nGDOF = 3*i-2;
			SuppressDOF (nGDOF);

			for (int j=nGDOF ;j<=m_nDOF;j++)				//Changing the structural load matrix according to B.C
			{ 
				if (j==nGDOF)
				{
					m_dSNF(j,1)=fXD;
				}

				else
				{
					m_dSNF(j,1)=m_dSNF(j,1)-m_dSSM(j,nGDOF)*fXD;
				}
			}
		}
		if (nYFC == 1)
		{
			int nGDOF = 3*i-1;
			SuppressDOF (nGDOF);
			for (int k=nGDOF ;k<=m_nDOF;k++)			//changing the structural load matrix according to B.C
			{ if (k==nGDOF)
			{
				m_dSNF(k,1)=fYD;
			}
			else
			{
				m_dSNF(k,1)=m_dSNF(k,1)-m_dSSM(k,nGDOF)*fYD;
			}
			}

		}
		if (nZFC == 1)
		{
			int nGDOF = 3*i;
			SuppressDOF (nGDOF);
			for (int n=nGDOF ;n<=m_nDOF;n++)            //changing the structural load matrix according to B.C
			{ if (n==nGDOF)
			m_dSNF(n,1)=fZD;
			else
				m_dSNF(n,1)=m_dSNF(n,1)-m_dSSM(n,nGDOF)*fZD;
			}
		}
	}

	// debug?
	if (m_nDebugLevel == 1)
	{
		MTB.PrintMatrixRowWise (m_dSSM, "Structural Stiffness (After BCs)",     
			m_FileOutput);
		MTB.PrintMatrixColumnWise (m_dSNF, "Structural Nodal Forces (After BCs)",   
			m_FileOutput);
	}
}

void CTruss::SuppressDOF (const int nEqn)                     
	// ---------------------------------------------------------------------------
	// Function: Imposes the essential boundary conditions
	// Input:    Equation number
	// Output:   none
	// ---------------------------------------------------------------------------
{

	for (int j=1; j <= m_nDOF; j++)
	{
		// zero out the row
		m_dSSM(nEqn, j) = 0.0;
		// zero out the column
		m_dSSM(j, nEqn) = 0.0;
	}
	// set diagonal to 1  
	m_dSSM(nEqn, nEqn) = 1.0;
}
void CTruss::Solve ()
	// ---------------------------------------------------------------------------
	// Function: Solves the system equations for the nodal displacements
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	double TOL = 1.0e-16;
	if(MTB.LDLTFactorization (m_dSSM, TOL) )       //Factorizing the stiffness matrix
	{
		
		if(MTB.LDLTSolve(m_dSSM,m_dSND,m_dSNF))     //solving for the structural nodal displacements
		{
			for ( int i=1; i <= m_nNodes; i++)
			{
				float m_fXDisp = static_cast<float>(m_dSND(3*i-2,1));
				float m_fYDisp = static_cast<float>(m_dSND(3*i-1,1));
				float m_fZDisp = static_cast<float>(m_dSND(3*i,1));
				m_NodalResponseData(i).SetDisplacements (m_fXDisp, m_fYDisp,m_fZDisp);
			}			
		}
	}
	
	//calculating nodal support reactions												CHANGED THIS


	for (int i=1;i<=m_nNodes;i++)                	
	{
		double xsr=0;
		double ysr=0;
		double zsr=0;
		float fx,fy,fz;
		m_NodalData(i).GetLoads( fx,fy,fz);

		if(fx==0)
		{
			for(int j=1;j<=m_nDOF;j++)
			{
				xsr+=m_dSSM(3*i-2,j)*m_dSND(j,1);
			}
		}
		else xsr=0;

		if(fy==0)
		{
			for(int j=1;j<=m_nDOF;j++)
			{
				ysr+=m_dSSM(3*i-1,j)*m_dSND(j,1);
			}
		}
		else ysr=0;

		if(fz==0)
		{
			for(int j=1;j<=m_nDOF;j++)
			{
				zsr+=m_dSSM(3*i,j)*m_dSND(j,1);
			}
		}
		else zsr=0;
		m_NodalResponseData(i).SetSupportReaction(xsr,ysr,zsr);

	}

}

void CTruss::Response ()
	// ---------------------------------------------------------------------------
	// Function: Computes the element response
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	int i, j;
	CVector<int> nE(6); // stores the element dof
	CMatrix<double> dT(2,6), dTT(6,2); // transformation matrices
	CMatrix<double> dLD(2,1), dND(6,1); // nodal displacements         
	float fX1, fX2, fY1, fY2, fL, fZ1, fZ2;
	float fA;
	double fE;
	int nSN, nEN;

	// initialize
	dT.Set (0.0);

	// loop thro' all elements
	for (i=1; i <= m_nElements; i++)
	{
		float TempDiff;
		double da;
		m_ElementData(i).GetData(nSN,nEN,fA,fE,da);

		//Calculating the temperature difference between element's start and end node. 

		TempDiff=m_NodalData(nEN).GetTemperature()-m_NodalData(nSN).GetTemperature();

		// form strain, stress and force in two steps
		m_NodalData(nSN).GetCoords (fX1, fY1, fZ1);
		m_NodalData(nEN).GetCoords (fX2, fY2, fZ2);
		fL = sqrt((fX2-fX1)*(fX2-fX1) +
			(fY2-fY1)*(fY2-fY1)+(fZ2-fZ1)*(fZ2-fZ1));
		// local-to-global transformation matrix

		double dl = double((fX2-fX1)/fL);
		double dm = double((fY2-fY1)/fL);
		double dn = double((fZ2-fZ1)/fL);
		dT(1,1) = dl; dT(1,2) = dm ; dT(1,3) = dn;
		dT(2,4) = dl; dT(2,5) = dm ; dT(2,6) = dn;

		// get element nodal displacements                       
		nE(1) = 3*nSN-2; nE(2) = 3*nSN-1;nE(3) = 3*nSN;
		nE(4) = 3*nEN-2; nE(5) = 3*nEN-1;nE(6) = 3*nEN;
		for (j=1; j <= 6; j++)
		{
			dND(j,1) = m_dSND(nE(j),1);
		}

		// form d'=T*d
		if(MTB. Multiply ( dT, dND,dLD))	                                 						
		{
			// strain
			float fStrain =float((dLD(2,1) - dLD(1,1))/fL)-(da*TempDiff);
			// stress
			float fStress = float(fE*fStrain);
			// force
			float fForce = float(fStress*fA);

			// update with the computed values
			m_ElementResponseData(i).SetData (fStrain,
				fStress, fForce);
		}
	}
}

void CTruss::CreateOutput ()
	// ---------------------------------------------------------------------------
	// Function: Creates the output file containing the results.
	//           Currently incomplete.
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	int i;
	float fX,fY,fZ;

	// print the problem size
	m_FileOutput << '\n';
	m_FileOutput << "PROBLEM SIZE" << '\n';
	m_FileOutput << "------------" << '\n';
	m_FileOutput << "   Number of nodes : " << m_nNodes << '\n';
	m_FileOutput << "Number of elements : " << m_nElements << '\n';

	// print the nodal coordinates
	m_FileOutput << '\n';
	m_FileOutput << "NODAL COORDINATES" << '\n';
	m_FileOutput << "-----------------" << '\n';
	m_FileOutput << "Node" << "   " << "X-Coordinate"
		<< "    " << "Y-Coordinate" << "    " << "Z-Coordinate"<< '\n';
	m_FileOutput << "----" << "   " << "------------"
		<< "    " << "------------"<< "    "  << "------------"<< '\n';
	for( i=1;i<=m_nNodes;i++)
	{
		m_NodalData(i).GetCoords( fX, fY,fZ);
		m_FileOutput<<  i    << "        " << fX
			<< "            " << fY << "            " <<fZ<< '\n';
	}

	// print the nodal fixities
	int nXFC,nYFC,nZFC;
	float fXD,fYD,fZD;
	m_FileOutput << '\n';
	m_FileOutput << "NODAL FIXITIES" << '\n';
	m_FileOutput << "--------------" << '\n';
	m_FileOutput << "Node" << "   " << "X-Fixity"
		<< "    " << "X-Disp" << "    " << "Y-Fixity"<< "    " << "Y-Disp"<< "    " << "Z-Fixity"<< "    " << "Z-Disp"<< '\n';
	m_FileOutput << "----" << "   " << "--------"
		<< "    " << "--------" << "    " << "--------" << "   " << "--------"<< "   " << "--------"<< "   " << "--------"<< '\n';
	for(i=1;i<=m_nNodes;i++)
	{
		m_NodalData(i).GetFixity( nXFC, nYFC,nZFC);
		m_NodalData(i).GetNodalSettlement(fXD,fYD,fZD);
		m_FileOutput<<  i    << "        " << nXFC<< "        " << fXD
			<< "            " << nYFC << "        " << fYD<< "        " << nZFC<< "        " << fZD<< '\n';
	}
	// print the nodal forces
	float fXL,fYL,fZL;
	m_FileOutput << '\n';
	m_FileOutput << "NODAL FORCES" << '\n';
	m_FileOutput << "------------" << '\n';
	m_FileOutput << "Node" << "   " << "X-Force" 
		<< "    " << "Y-Force"<< "    " << "Z-Force"  << "    " << "Temperature"  << '\n';
	m_FileOutput << "----" << "   " << "-------" 
		<< "    " << "-------"<< "   " << "-------" << "   " << "-------" << '\n';
	for(i=1;i<=m_nNodes;i++)
	{
		m_NodalData(i).GetLoads( fXL, fYL,fZL);  
		m_FileOutput<<  i    << "        " << fXL
			<< "            " << fYL << "        " << fZL<< "        " <<m_NodalData(i).GetTemperature()<< '\n';
	}

	//Element Data
	int nS,nE;
	float fA;
	double dE,dCTE;
	m_FileOutput << '\n';
	m_FileOutput << "ELEMENT DATA" << '\n';
	m_FileOutput << "--------------" << '\n';
	m_FileOutput << "ELEMENT" << "   " << "SN"
		<< "       " << "  EN" << "      " << "  Area"<< "     " << "  Elasticity"<< "      " << "  CTE"<<'\n';
	m_FileOutput << "----" << "   " << "--------"
		<< "    " << "--------" << "    " << "--------" << "   " << "--------"<< "   " << "--------"<<'\n';
	for(i=1;i<=m_nElements;i++)
	{
		m_ElementData(i).GetData( nS,nE,fA,dE,dCTE);
		m_FileOutput<<  i    << "        " << nS<< "        " << nE
			<< "            " << fA << "        " << dE<< "        " << dCTE<< '\n';
	}

	// print the nodal displacements
	m_FileOutput << '\n';
	float fXDisp, fYDisp,fZDisp;
	m_FileOutput << "NODAL DISPLACEMENTS" << '\n';
	m_FileOutput << "-------------------" << '\n';
	m_FileOutput << "Node" << "   " << "X-Displacement"
		<< "    " << "Y-Displacement"<< "    " << "Z-Displacement" << '\n';
	m_FileOutput << "----" << "   " << "--------------"
		<< "    " << "--------------"<< "   " << "--------------" << '\n';
	for (i=1; i <= m_nNodes; i++)
	{
		m_NodalResponseData(i).GetDisplacements (fXDisp, fYDisp, fZDisp);
		if(fXDisp!=0||fYDisp!=0||fZDisp!=0)
		{
			m_FileOutput << std::setw(4) << i << "   "
				<< std::setw(14) << fXDisp 
				<< "    " << std::setw(14) << fYDisp<< "    " << std::setw(14) << fZDisp << '\n';
		}
	}
	// print the element strain, stress and force
	float fStrain,fStress,fForce;
	m_FileOutput << '\n';
	m_FileOutput << "ELEMENT RESPONSE" << '\n';
	m_FileOutput << "----------------" << '\n';
	m_FileOutput << "  Element" << "    " << " Strain" 
		<< "       " << " Stress" << "        " << "   Force" << '\n';
	m_FileOutput << "-------" << "        " << "------" << "       " 
		<< "------" << "       " << "-----" << '\n';
	for (i=1; i <= m_nElements; i++)
	{
		m_ElementResponseData(i).GetData (fStrain, fStress,fForce);
		m_FileOutput << std::setw(4) << i << "   "
			<< std::setw(10) << fStrain 
			<< "     " << std::setw(13) << fStress<< std::setw(13)<<"  " << fForce << '\n';
	}
	
	//Print nodal reactions																		CHANGED THIS					
	double xsr,ysr,zsr;
	m_FileOutput<<'\n';
	m_FileOutput<<"                              Nodal Support Reactions"<<'\n';
	m_FileOutput<<"                              -----------------------"<<'\n';
	m_FileOutput<<"Node"<<"   "<<"X Reaction"<<"   "<<"Y Reaction"<<"   "<<"Z Reaction"<<'\n';
	m_FileOutput<<"----"<<"   "<<"----------"<<"   "<<"----------"<<"   "<<"----------"<<'\n';
	for(int i=1;i<=m_nNodes;i++)
	{
		m_NodalResponseData(i).GetSupportReaction(xsr,ysr,zsr);
		m_FileOutput<<i<<"    "<< std::setw(10)<<xsr <<"    "<< std::setw(10)<<ysr<<"    "<< std::setw(10)<<zsr<<'\n';

	}
	m_FileOutput<<"\n\nAbsolute error:"<<m_fAE<<'\n';
	m_FileOutput<<"\n\nRelative Error:"<<m_fRE<<'\n';
	

}


//Calculating the Absolute and relative error norms

void CTruss::GetErrors()
{

	CMatrix<double>dTemp1;
	CMatrix<double>dTemp2;
	double sum1=0;
	double sum2=0;
	double N;
	
	if(MTB.Multiply(m_dSSM,m_dSND,dTemp1))
	{
		if(MTB.Subtract(dTemp1,m_dSNF,dTemp2))
		{for(int i=1;i<=m_nNodes;i++)
		{
			sum1+=dTemp2(i,1)*dTemp2(i,1);

		}
		m_fAE=float(sqrt(sum1));
		for(int i=1;i<=m_nNodes;i++)
		{
			sum2+=m_dSNF(i,1)*m_dSNF(i,1);
		}
		N=sqrt(sum2);
		m_fRE=m_fAE/N;

		}
	}
    
}

void CTruss::TerminateProgram ()
	// ---------------------------------------------------------------------------
	// Function: Closes input and output files
	// Input:    none
	// Output:   none
	// ---------------------------------------------------------------------------
{
	// close the input and output files
	m_FileInput.close ();
	m_FileOutput.close ();

	std::cout << "\nExecution completed successfully." << std::endl;
}

void CTruss::ErrorHandler (ERRORCODE nCode) const
	// ---------------------------------------------------------------------------
	// Function: Displays the error message on the error stream
	// Input:    Error #
	// Output:   none
	// ---------------------------------------------------------------------------
{
	std::cerr << '\n';
          if (nCode == INVALIDNUMNODES)                                 // invalid number of nodes
		std::cerr << "Number of nodes must be >= 2.";
	else if (nCode == INVALIDNUMELEMENTS)                         // invalid number of elements
		std::cerr << "Number of elements must be >= 1.";
	else if (nCode == INVALIDDEBUGCODE)                           // invalid debug level
		std::cerr << "Debug level must be 0 or 1.";
	else if (nCode == INVALIDNODENUM)                             // invalid node number
		std::cerr << "Invalid node number";
	else if (nCode == INVALIDELEMENTNUM)                          // invalid element number
		std::cerr << "Invalid element number";
	else if (nCode == INVALIDCSAREA)                              // invalid x/s area
		std::cerr << "Area must be positive.";
	else if (nCode == INVALIDYM)                                  // invalid E
		std::cerr << "Modulus of elasticity must be positive.";
	else if (nCode == UNSTABLETRUSS)                              // unstable structure
		std::cerr << "Unstable truss.";
	else if (nCode == INVALIDINPUT)                               // invalid input
		std::cerr << "Input file need *end as last line.";
	else if (nCode == INVALIDNODALFIXITY)                         // invalid fixity code
		std::cerr << "Nodal fixity code must be 0 or 1.";
	else if (nCode == MISSINGEND)                                 // missing end statement
		std::cerr << "Missing *END statement in input file.";
	else if (nCode == INVALIDCOMMANDLINE)                         // need 1 or 3 command line arguments
		std::cerr << "Invalid number of command line arguments.";
	else if(nCode==INDETERMINANT)																//CHANGED THIS
		std::cerr<<"Truss Indeterminant";
	if (nCode != UNSTABLETRUSS && nCode != INVALIDCOMMANDLINE)
		std::cerr << '\n' << "Error in input file line : " 
		<< m_nLineNumber;

	std::cerr << std::endl;

	// fatal error, exit the program
	exit (1);
}