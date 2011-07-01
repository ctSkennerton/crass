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
#include <map>


typedef std::pair<std::string, std::string> stringPair;
double similarityScore(char a, char b);
double findMax(double array[], int length);
stringPair smithWaterman(std::string seqA, std::string seqB);

#endif
