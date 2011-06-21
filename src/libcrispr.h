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
#include "CrisprNode.h"
#include "NodeManager.h"
#include "WuManber.h"
#include "libbitap.h"
#include "bm.h"
#include "kseq.h"
#include "ReadHolder.h"
#include "SeqUtils.h"


typedef std::map<std::string, CrisprNode *> lookupTable;

typedef std::vector<ReadHolder *> ReadList;
typedef std::vector<ReadHolder *>::iterator ReadListIterator;

// direct repeat as a string and a list of the read objects that contain that direct repeat
typedef std::map<std::string, ReadList *> ReadMap;
typedef std::map<std::string, ReadList *>::iterator ReadMapIterator;

// The following two classes are simple data storages objects
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
        int  DR_Length;
        int  DR_StartPos;                   // the start of the 'right' dr in the read
        int  DR_EndPos;                     // the end of the whole direct repeat
        int  DR_MatchStartPos;              // the start of the 'left' dr in the read
        int  DR_MatchEndPos;
        int  DR_NumMismatches;              // difference between either string
};


class ReadMatch {
    public:
        // constructor / destructor
        ReadMatch ();
        ~ReadMatch () {}
        
        // members
        const char * RM_SubstrStart;        // the read sequence
        const char * RM_SubstrEnd;          // the read sequence minus the start to the pattern
        int  RM_StartPos;                   // the start position
        int  RM_EndPos;                     // the end pos
        int  RM_MatchStartPos;              // the start position
        int  RM_MatchEndPos;                // the end pos
        int  RM_NumMismatches;
        int  RM_NumInsertions;
        int  RM_NumDeletions;
        int  RM_NumSubstitutions;
};

//**************************************
// search functions
//**************************************

// boyer moore functions

float bm_search_fastq_file(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &spacerLookup, lookupTable &kmerLookup, ReadMap * mReads);

bool inline check_mismatch( int &temp_mismatch, const options &opts);

// bitap functions
bool direct_repeat_init( DirectRepeat &dr_match, ReadMatch &match_info, int &temp_mismatch, const options &opts );

bool direct_repeat_extend ( DirectRepeat &dr_match, ReadMatch &match_info, int &temp_mismatch, const options &opts);

bool update_word_bitap(ReadMatch &match_info, int &search_begin, int &start, const options &opts, DirectRepeat &dr, std::string &subject_word, int &temp_mismatch);

float bitap_search_fastq_file(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &spacerLookup, lookupTable &kmerLookup, ReadMap *mReads);

// multimap functions (wuMander)
void read_for_multimatch(const char *input_fastq, const options &opts, std::vector<std::string> &patterns, lookupTable &directRepeatLookup, ReadMap *mReads );

//**************************************
// kmer operators
//**************************************

void cutLeftKmer( std::string &read, int &start, int &end, lookupTable &inputLookup, const options &opts);

void cutRightKmer( std::string &read, int &start, int &end, lookupTable &inputLookup, const options &opts);

void cutSpacerKmers( std::string &read, lookupTable &inputLookup, const options &opts);

bool cut_direct_repeat_sequence(DirectRepeat &dr_match, const options &opts, string &read);

bool check_dr_and_spacer_length(const options &opts, DirectRepeat &dr_match);


//**************************************
// Read Holder
//**************************************

void DRLowLexi(std::string matchedRead, DirectRepeat * dr_match,  ReadHolder * tmp_holder);

void addReadHolder(ReadMap * mReads, ReadHolder * tmp_holder, DirectRepeat * dr_match, std::string read_header, std::string read);

//**************************************
// Trimming
//**************************************

void drTrim();



//**************************************
// lookup table shite
//**************************************

void addToLookup(const string &dr ,lookupTable &patternsLookup);

//**************************************
// system
//**************************************

gzFile getFileHandle(const char * inputFile);

//**************************************
// STL extensions
//**************************************

bool inline keyExists(lookupTable &patterns_hash, std::string &direct_repeat);

bool inline keyExists(ReadMap * mReads, std::string &direct_repeat);

void map_to_vector(lookupTable &patterns_hash, std::vector<std::string> &patterns);

#endif //libcrispr_h