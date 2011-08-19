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
#include "crass_defines.h"
#include "WuManber.h"
#include "libbitap.h"
#include "bm.h"
#include "kseq.h"
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "StringCheck.h"
#include "Genome.h"

typedef std::map<std::string, bool> lookupTable;

typedef std::vector<ReadHolder *> ReadList;
typedef std::vector<ReadHolder *>::iterator ReadListIterator;

// direct repeat as a string and a list of the read objects that contain that direct repeat
typedef std::map<StringToken, ReadList *> ReadMap;
typedef std::map<StringToken, ReadList *>::iterator ReadMapIterator;

enum READ_TYPE {
    LONG_READ,
    SHORT_READ,
};

enum side{rightSide, leftSide};

// The following  class is simple data storages objects
// they are stupidly public for that reason

class DirectRepeat {
    //-----
    // 
    //
    public:
        // constructor / destructor
        DirectRepeat();
        ~DirectRepeat() {}
        
        void reset(void);

        
    // members
        std::string DR_Sequence;
        std::string DR_MatchSequence;
        std::string DR_Spacer;
        std::vector<int> DR_StartStopList;      // a vector containing the starting positions of all the DR
        int  DR_Length;
        int  DR_StartPos;                   // the start of the 'right' dr in the read
        int  DR_EndPos;                     // the end of the whole direct repeat
        int  DR_MatchStartPos;              // the start of the 'left' dr in the read
        int  DR_MatchEndPos;
        int  DR_NumMismatches;              // difference between either string
};



//**************************************
// search functions
//**************************************

READ_TYPE decideWhichSearch(const char *input_fastq);

float crtSearchFastqFile(const char *input_fastq, const options &opts, ReadMap * mReads, StringCheck * mStringCheck);

float bmSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck);

void findSingletons(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap *mReads, StringCheck * mStringCheck);

bool partialStarting (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq, int * f_start, int * f_end);

bool partialEnding (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq, int * f_start, int * f_end);
//**************************************
// kmer operators
//**************************************

void cutLeftKmer( std::string &read, int &start, int &end, lookupTable &inputLookup, const options &opts);

void cutRightKmer( std::string &read, int &start, int &end, lookupTable &inputLookup, const options &opts);

void cutSpacerKmers( std::string &read, lookupTable &inputLookup, const options &opts);

bool cutDirectRepeatSequence(DirectRepeat &dr_match, const options &opts, string &read);

bool checkDRAndSpacerLength(const options &opts, DirectRepeat &dr_match);

bool isLowComplexity(DirectRepeat &dr_match);

bool isSpacerAndDirectRepeatSimilar(DirectRepeat &dr_match);

//int getActualRepeatLength(std::vector<int> &candidateCRISPR, std::string &read, int searchWindowLength, int minSpacerLength);

//**************************************
// Read Holder
//**************************************

std::string DRLowLexi(std::string& matchedRead,  ReadHolder * tmp_holder);

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmp_holder, std::string& read_header, std::string& read);

//**************************************
// lookup table shite
//**************************************

void addToLookup(const string &dr ,lookupTable &patternsLookup);



//**************************************
// STL extensions
//**************************************

bool inline keyExists(lookupTable &patterns_hash, std::string &direct_repeat);

void map2Vector(lookupTable &patterns_hash, std::vector<std::string> &patterns);

//*************************************
// Old genome finder stuff
//*************************************
void printGenomeCrispr(std::vector<Crispr*>& CRISPRVector, options& opts, std::string& mHeader, std::string& mSequence, bool repeatsFound);

bool  findRepeats(void);

void  trim(Crispr * candidateCRISPR, int minRepeatLength);

void  checkFlank(side sT, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence,int sequenceLength, std::string& read);

int scan(side sT, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double confidence, int sequenceLength, std::string& sequence);

void scanRight(Crispr * candidateCRISPR, std::string& pattern, int minSpacerLength, int scanRange, int sequenceLength, std::string& sequence);

bool  patternMatches(std::string& pattern1, std::string& pattern2, double confidence);

int   min (int * array);

int   getHammingDistance(std::string& seq1, std::string& seq2);

#endif //libcrispr_h
