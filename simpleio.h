/*********************************************
Utility Library Function
Copyright(c) 2000, S. D. Rajan
All rights reserved

Introduction to Structural Analysis and Design
*********************************************/
#ifndef __RAJAN_SIMPLEIO_H__
#define __RAJAN_SIMPLEIO_H__

#pragma once

#include <string>
using std::string;

#include <iostream>
using std::cin;
using std::cout;

template <typename T>
void GetInteractive (const string& szPrompt, T& Value)
{
	int nState;
	string szTemp;

	do {
		cout << szPrompt;
		cin >> Value;
		nState = cin.fail(); // state is non-zero if cin fails
		if (nState) {
			cin.clear (0);
			cin >> szTemp;	// since cin failed store the invalid
							// input in szTemp
		}
	} while (nState);	// iterate until the input is valid.
						// state is zero for a valid input
}

#endif