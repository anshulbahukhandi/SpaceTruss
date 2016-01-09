/*********************************************
Planar Truss Analysis Program
Copyright(c) 2000-13, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
Object-Oriented Numerical Analysis
*********************************************/
#include "truss.h"

int main (int argc, char* argv[])
{
	CTruss MyTruss; // the one and only truss

	// show program banner
	MyTruss.Banner (std::cout);

	// open input and output files
	MyTruss.PrepareIO (argc, argv);
	
	// read truss model
	MyTruss.ReadTrussModel ();

    // analyze
    MyTruss.Analyze ();

	// close input and output files
	MyTruss.TerminateProgram ();

	return 0;
}