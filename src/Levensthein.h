//
//  Levensthein.h
//  MCD
//
//  Created by Connor Skennerton on 15/05/11.
//  Copyright 2011 The Faculty of EAIT. All rights reserved.
//
#ifndef __LEV_H
#define __LEV_H
#include <string>
#include <vector>
using namespace std;

typedef vector< vector<int> > Tmatrix; 


int Levensthein_distance( string source,  string target);
#endif
