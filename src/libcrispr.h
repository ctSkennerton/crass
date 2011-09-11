// File: libcrispr.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Header file for the "crispr toolbox"
//
// --------------------------------------------------------------------
//  Copyright  2011 Michael Imelfort and Connor Skennerton
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//

#ifndef libcrispr_h
#define libcrispr_h

// system includes
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h> 
#include <map>

// local includes
#include "crassDefines.h"
#include "WuManber.h"
#include "bm.h"
#include "kseq.h"
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "StringCheck.h"

typedef std::map<std::string, bool> lookupTable;

typedef std::vector<ReadHolder *> ReadList;
typedef std::vector<ReadHolder *>::iterator ReadListIterator;

// direct repeat as a string and a list of the read objects that contain that direct repeat
typedef std::map<StringToken, ReadList *> ReadMap;
typedef std::map<StringToken, ReadList *>::iterator ReadMapIterator;



enum READ_TYPE {
    LONG_READ,
    SHORT_READ
};

enum side{rightSide, leftSide};


//**************************************
// search functions
//**************************************

READ_TYPE decideWhichSearch(const char *inputFastq, float * aveReadLength);

void longReadSearch(const char *input_fastq, const options &opts, ReadMap * mReads, StringCheck * mStringCheck);

void shortReadSearch(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck);

void findSingletons(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap *mReads, StringCheck * mStringCheck);


//**************************************
// Read Holder
//**************************************

std::string DRLowLexi( ReadHolder * tmp_holder, std::string& read);

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmp_holder);

//**************************************
// lookup table shite
//**************************************

void addToLookup(const string &dr ,lookupTable &patternsLookup);

float getStringSimilarity(std::string& s1, std::string& s2);

//**************************************
// STL extensions
//**************************************

bool inline keyExists(lookupTable &patterns_hash, std::string &direct_repeat);

void map2Vector(lookupTable &patterns_hash, std::vector<std::string> &patterns);

//*************************************
// Ported CRT code
//*************************************

unsigned int extendPreRepeat(ReadHolder* tmp_holder, int searchWindowLength, int minSpacerLength);

bool qcFoundRepeats(ReadHolder * tmp_holder);

bool isRepeatLowComplexity(std::string& repeat);

int scanRight(ReadHolder * tmp_holder, std::string& pattern, unsigned int minSpacerLength, unsigned int scanRange);

bool  patternMatches(std::string& pattern1, std::string& pattern2, double confidence);


#endif //libcrispr_h
