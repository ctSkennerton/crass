//
//  SmithWaterman.h
//  crass
//
//  Created by Connor Skennerton on 23/06/11.
//  Copyright 2011 Australian Centre for Ecogenomics. All rights reserved.
//

#ifndef __SMITH_WATERMAN_H
#define __SMITH_WATERMAN_H
#include <string>

double similarityScore(char a, char b);
double findMax(double array[], int length);
int smithWaterman(std::string seqA, std::string seqB );

#endif
