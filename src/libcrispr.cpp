// File: libcrispr.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// This file wraps all the crispr associated functions we'll need.
// The "crispr toolbox"
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
// system includes
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <zlib.h>  
#include <fstream>
#include <fcntl.h>
#include <stdlib.h>


// local includes
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "CrisprNode.h"
#include "NodeManager.h"
#include "WuManber.h"
#include "libbitap.h"
#include "bm.h"
#include "SeqUtils.h"
#include "Levensthein.h"

//**************************************
// DirectRepeat
//**************************************
DirectRepeat::DirectRepeat() {
    //-----
    // Constructor!
    //
    DR_Sequence = "";
    DR_MatchSequence = "";
    DR_Spacer = "";
    DR_Length = 0;
    DR_MatchStartPos = 0;
    DR_MatchEndPos = 0;
    DR_StartPos = 0;
    DR_EndPos = 0;
    DR_NumMismatches = 0;
}

void DirectRepeat::reset(void)
{
    //-----
    // DER
    //
    DR_Sequence = "";
    DR_MatchSequence = "";
    DR_Spacer = "";
    DR_Length = 0;
    DR_MatchStartPos = 0;
    DR_MatchEndPos = 0;
    DR_StartPos = 0;
    DR_EndPos = 0;
    DR_NumMismatches = 0;
}
/*
**************************************
 RepeatMatch
**************************************
ReadMatch::ReadMatch() {
    //-----
    // Constructor!
    //
    RM_SubstrStart = NULL;
    RM_SubstrEnd = NULL;
    RM_StartPos = 0;
    RM_EndPos = 0;
    RM_MatchStartPos = 0;
    RM_MatchEndPos = 0;
    RM_NumMismatches = 0;
    RM_NumInsertions = 0;
    RM_NumDeletions = 0;
    RM_NumSubstitutions = 0;
}
*/
/* 
declare the type of file handler and the read() function
as described here:
http://lh3lh3.users.sourceforge.net/parsefastq.shtml

THIS JUST DEFINES A BUNCH OF **templated** structs

*/
KSEQ_INIT(gzFile, gzread);  

//**************************************
// search functions
//**************************************

// boyer moore functions
float bmSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads)
{
    gzFile fp = getFileHandle(input_fastq);
    kseq_t *seq;
    int l, match_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    DirectRepeat dr_match;
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        std::string read = seq->seq.s;
        std::string read_header = seq->name.s;
        int seq_length = read.length() - 1;
        int search_end = seq_length - opts.lowDRsize;
        
        logInfo("read counter: "<<match_counter, 8);
        
        total_base += seq_length;
        
        bool match_found = false;
        // create the read holder
        ReadHolder * tmp_holder = new ReadHolder();
        
        // boyer-moore search
        for (int start = 0; start < search_end; start++)
        {
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            logInfo("search begin: "<<search_begin<<" search end: "<<search_end, 8);
            
            if (search_begin >= search_end ) break;
            
            
            std::string query_word = read.substr(start, opts.lowDRsize);
            std::string subject_word = read.substr(search_begin);
            
            logInfo("query: "<<query_word<<" subject: "<<subject_word, 10);
            
            
            int MatchStartPos = PatternMatcher::bmpSearch( subject_word, query_word );
            logInfo("bm return: "<<MatchStartPos, 8);
            
            if (MatchStartPos > -1) 
            {
                int EndPos = MatchStartPos + search_begin + opts.lowDRsize;
                int MatchEndPos = start + opts.lowDRsize;
                
                dr_match.DR_StartPos = MatchStartPos + search_begin;
                //dr_match.DR_StartList.push_back(MatchStartPos + search_begin);
                
                dr_match.DR_MatchStartPos = start;
                //dr_match.DR_StartList.push_back(start);
                
                dr_match.DR_MatchEndPos = MatchEndPos;
                dr_match.DR_EndPos = EndPos;
                
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                ++EndPos;
                ++MatchEndPos;
                if (EndPos <= seq_length) 
                {
                    // read through the subsuquent bases untill they don't match
                    logInfo("Read: "<<seq->name.s<<" Len: "<<read.length(), 10);
                    while (read.at(MatchEndPos) == read.at(EndPos)) 
                    {
                        logInfo("Match end pos: "<<MatchEndPos<<" end pos: "<<EndPos, 10);
                        logInfo(read.at(MatchEndPos) << " : " << MatchEndPos << " == " << read.at(EndPos) << " : " << EndPos, 10);

                        
                        dr_match.DR_MatchEndPos = MatchEndPos;
                        dr_match.DR_EndPos = EndPos;
                        ++EndPos;
                        ++MatchEndPos;
                        if (EndPos > seq_length) break;
                    }
                    
                    if (cutDirectRepeatSequence(dr_match, opts, read))
                    {
                        patterns_hash[dr_match.DR_MatchSequence] = true;

                        
                        tmp_holder->RH_StartStops.push_back( dr_match.DR_MatchStartPos);
                        tmp_holder->RH_StartStops.push_back( dr_match.DR_MatchEndPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                        match_found = true;
                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
                        continue;
                        /*
                        // check if there is alot of read still to go ( long read )
                        // if the end pos is less than or equal to the end of the read minus 2 * the 
                        // minimun direct repeat size and the minimum spacer size
                        if (dr_match.DR_EndPos <= (search_end - opts.lowDRsize - opts.lowSpacerSize))
                        {
                            
                            // if yes then use the full sequence to find more instances of this string
                            // cut a substring minusing the first and last two bases ( in case of a overextend )
                            std::string new_search_string = dr_match.DR_Sequence.substr(2, dr_match.DR_Sequence.length() - 4);
                            
                            std::vector<int> start_list;

                            // get all of the positions of this substring
                            PatternMatcher::bmpMultiSearch(read.substr(dr_match.DR_EndPos), new_search_string, start_list);
                            
                            
                            std::vector<int>::iterator vec_iter = start_list.begin();
                            
                            while (vec_iter != start_list.end()) 
                            {
                               *vec_iter = *vec_iter + dr_match.DR_EndPos;
                                ++vec_iter;
                            }
                            vec_iter = start_list.begin();
                            start_list.insert(vec_iter, (dr_match.DR_StartPos + 2));

                            vec_iter = start_list.begin();
                            start_list.insert(vec_iter, (dr_match.DR_MatchStartPos + 2));

                            vec_iter = start_list.begin();

                            
                            while (vec_iter != start_list.end()) 
                            {
                                logInfo("start pos of multi search: "<<*vec_iter, 10);
                                ++vec_iter;
                            }
                            // update the start positions of the direct repeats 
                            int real_length = getActualRepeatLength(start_list, read, new_search_string.length(), opts.lowSpacerSize);
                            
                            logInfo("DR length: "<<real_length, 10);
                            
                            vec_iter = start_list.begin();
                            while (vec_iter != start_list.end()) 
                            {
                                dr_match.DR_StartStopList.push_back(*vec_iter);
                                dr_match.DR_StartStopList.push_back(real_length);
                                ++vec_iter;
                            }
                            addReadHolder(mReads, tmp_holder, &dr_match, read_header, read);
                            break;
                            
                        }
                        else
                        {
                            
                            dr_match.DR_StartStopList.push_back(dr_match.DR_MatchStartPos);
                            dr_match.DR_StartStopList.push_back(dr_match.DR_MatchEndPos);
                            dr_match.DR_StartStopList.push_back(dr_match.DR_StartPos);
                            dr_match.DR_StartStopList.push_back(dr_match.DR_EndPos);
                            // if no then add the read holder to the map
                            addReadHolder(mReads, tmp_holder, &dr_match, read_header, read);
                        
                            start = dr_match.DR_StartPos - 1;
                            dr_match.reset();
                            continue;
                        }
                         */
                    } 
                    else 
                    {                        
                        // increment the for loop so that it begins at the start of the 
                        // matched kmer/direct repeat
                        // minus 1 cause it will be incremented again at the top of the for loop
                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
                        continue;
                    }
                }
            }
        }
        if (match_found)
        {
            readsFound[read_header] = true;
            addReadHolder(mReads, tmp_holder, read_header, read);
        }

        match_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<input_fastq, 1);
    return total_base / match_counter;
}


float bitapSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap *mReads) 
{
    
    gzFile fp = getFileHandle(input_fastq);
    kseq_t *seq;
    int l, match_counter = 0;
    unsigned int total_base = 0;
    
    // initialize seq
    seq = kseq_init(fp);
    
    DirectRepeat dr_match;
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        int seq_length = (int) (strlen(seq->seq.s));
        int search_end = seq_length - opts.lowDRsize;
        bool match_found = false;

        total_base += seq_length;
        
        dr_match.reset();
        
        std::string read = seq->seq.s;

        
        bitapType b;
        
        ReadHolder * tmp_holder = new ReadHolder;

        for (int start = 0; start < search_end; start++)
        {
            // don't search too far into the read if we don't need to
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            logInfo("search begin: "<<search_begin<<" search end: "<<search_end, 10);
            
            if (search_begin >= search_end )
            {
                break;
            }
            // make sure that we are still under the number of mismatches
            if (dr_match.DR_NumMismatches > opts.max_mismatches)
            {
                break;
            }
            
            //stupid std::string concatenation stuff cause of the requirements of libbitap!!
            std::string substring_word = read.substr(start, opts.lowDRsize);
            std::string query_word = ".*" + substring_word + ".*";
            
            
            std::string subject_word = read.substr(search_begin, (read.length() - search_begin));
            //const char * 
            
            logInfo("query: "<<query_word<<" subject: "<<subject_word, 10);
            
            NewBitap(&b, query_word.c_str());
            
            
            const char * substrEnd = NULL;
            const char * substrStart = NULL;
            int numMismatches = 0;
            
            
            if (NULL != (substrEnd = FindWithBitap (&b, subject_word.c_str(), (int)subject_word.length(), opts.max_mismatches, &numMismatches, &substrStart))) 
            {
                dr_match.DR_MatchStartPos = start;
                dr_match.DR_StartPos = (int) (substrStart - subject_word.c_str()) + search_begin;

                
                dr_match.DR_MatchEndPos = start + opts.lowDRsize;
                substrEnd = substrStart + opts.lowDRsize;
                dr_match.DR_EndPos = dr_match.DR_StartPos + opts.lowDRsize ;
                
                            
                ++dr_match.DR_MatchEndPos;
                ++dr_match.DR_EndPos;
                
                while (numMismatches <= opts.max_mismatches && dr_match.DR_EndPos < seq_length) 
                {
                    char dr_A = read.at(dr_match.DR_MatchEndPos);
                    char dr_B = read.at(dr_match.DR_EndPos);
                    
                    logInfo("Match end pos: "<<dr_match.DR_MatchEndPos<<" end pos: "<<dr_match.DR_EndPos, 10);
                    logInfo(read.at(dr_match.DR_MatchEndPos) << " : " << dr_match.DR_MatchEndPos << " == " << read.at(dr_match.DR_EndPos) << " : " << dr_match.DR_EndPos, 10);
                    
                    if (dr_A != dr_B)
                    {
                        numMismatches++;
                        if (numMismatches > opts.max_mismatches)
                        {
                            // since the number of mismatches would have been one higher from the above loop
                            numMismatches--;
                            break;
                        }
                    }
                    
                    ++dr_match.DR_MatchEndPos;
                    ++dr_match.DR_EndPos;
                }
            
                dr_match.DR_NumMismatches = numMismatches; 
                
                if (cutDirectRepeatSequence(dr_match, opts, read))
                {

                    patterns_hash[dr_match.DR_MatchSequence] = true;
                    
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_MatchStartPos);
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_MatchEndPos);
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                    match_found = true;
                    start = dr_match.DR_StartPos - 1;
                    dr_match.reset();
                    continue;
                } 
                else 
                {
                    start = dr_match.DR_StartPos - 1;
                    dr_match.reset();
                    continue;
                }
            }
            DeleteBitap (&b);
        }
        
        if (match_found)
        {
            readsFound[seq->name.s] = true;
            addReadHolder(mReads, tmp_holder, seq->name.s, read);
        }
        ++match_counter;

    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<input_fastq, 1);

    return total_base / match_counter;
}

void scanForMultiMatches(const char *input_fastq, const options &opts, lookupTable &patterns_hsah, lookupTable &readsFound, ReadMap * mReads )
{
    std::vector<std::string> patterns;
    
    map2Vector(patterns_hsah, patterns);
    
    if (patterns.size() == 0)
    {
        logError("No patterns in vector for multimatch");
    }
    
    gzFile fp = getFileHandle(input_fastq);
    kseq_t *seq;
    //ReadMatch match_info;
    
    seq = kseq_init(fp);
    
    WuManber search;
    search.Initialize(patterns);
    DirectRepeat dr_match;
    
    int l;
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        dr_match.reset();
        
        //initialize with an impossible number
        int endPos = -1;
        
        std::string read = seq->seq.s;
        
        dr_match.DR_Sequence = search.Search(strlen(seq->seq.s), seq->seq.s, patterns, endPos);
        
        dr_match.DR_StartPos = endPos;
        if (endPos != -1)
        {
            std::string header = seq->name.s;
            if (!(keyExists(readsFound, header)))
            {
                dr_match.DR_EndPos = endPos + dr_match.DR_Sequence.length();
                
                // TODO: change this to some smarter logic
                // really silly way or breaking it up!!
                int thirds = read.length()/3;
                
                
                // create the read holder
                ReadHolder * tmp_holder = new ReadHolder;
                
                if (dr_match.DR_EndPos <= thirds)
                {
                    // first third
                    
                    // the match we already have would come first
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_StartPos);
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_EndPos);
                    partialEnding(dr_match, tmp_holder, read);
                    //findLostSouls(dr_match, tmp_holder, read, -1 );
                }
                else if (dr_match.DR_EndPos >= thirds*2)
                {
                    // last third
                    
                    partialStarting(dr_match, tmp_holder, read);
                    
                    // the match we already have would come last
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_StartPos);
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_EndPos);
                }
                else
                {
                    // middle third
                    if (partialStarting(dr_match, tmp_holder, read))
                    {
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                    }
                    else
                    {
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                        
                        partialEnding(dr_match, tmp_holder, read);
                    }
                }
                addReadHolder(mReads, tmp_holder, seq->name.s, read);
            }
        }
    }
    logInfo("finished multi pattern matcher", 1);
}

//**************************************
// kmer operators
//**************************************

bool cutDirectRepeatSequence(DirectRepeat &dr_match, const options &opts, string &read)
{
    dr_match.DR_Length = dr_match.DR_EndPos - dr_match.DR_StartPos;
    dr_match.DR_Sequence = read.substr(dr_match.DR_StartPos, (dr_match.DR_EndPos - dr_match.DR_StartPos + 1));
    dr_match.DR_MatchSequence = read.substr(dr_match.DR_MatchStartPos, (dr_match.DR_MatchEndPos - dr_match.DR_MatchStartPos + 1));
    dr_match.DR_Spacer = read.substr(dr_match.DR_MatchEndPos, (dr_match.DR_StartPos - dr_match.DR_MatchEndPos));

    if (!(checkDRAndSpacerLength(opts, dr_match)) || isLowComplexity(dr_match) || isSpacerAndDirectRepeatSimilar(dr_match))
    {
        return false;
    }
    
    return true;
}

bool checkDRAndSpacerLength(const options &opts, DirectRepeat &dr_match)
{
    logInfo("dr end pos: "<<dr_match.DR_EndPos<<" dr start pos: "<<dr_match.DR_StartPos, 10);
    
    int spacer_length = (int)dr_match.DR_Spacer.length();
    logInfo("DR len: "<<dr_match.DR_Length<<" SP length: "<<spacer_length, 10);
    // check if the direct repeat is in the right size range
    if ((dr_match.DR_Length < opts.lowDRsize) or (dr_match.DR_Length > opts.highDRsize)) 
    {
        return false;
    }
    
    // check if the spacer is in the right size range
    if ((spacer_length < opts.lowSpacerSize) or (spacer_length > opts.highSpacerSize))
    { 
        return false;
    }

    return true; 
}

bool isLowComplexity(DirectRepeat &dr_match)
{
    int cCount = 0;
    int gCount = 0;
    int aCount = 0;
    int tCount = 0;
    
    float aPercent;
    float cPercent;
    float gPercetn;
    float tPercent;

    std::string::iterator dr_iter = dr_match.DR_Sequence.begin();
    //int i = dr_match.DR_StartPos;
    while (dr_iter != dr_match.DR_Sequence.end()) 
    {
        switch (*dr_iter) 
        {
            case 'c':
            case 'C':
                cCount++; break;
            case 't': 
            case 'T':
                tCount++; break;
            case 'a':
            case 'A':
                aCount++; break;
            case 'g':
            case 'G':
                gCount++; break;
            default: break;
        }
        dr_iter++;
    }
    aPercent = aCount/dr_match.DR_Length;
    tPercent = tCount/dr_match.DR_Length;
    gPercetn = gCount/dr_match.DR_Length;
    cPercent = cCount/dr_match.DR_Length;
    
    if (aPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || tPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || gPercetn > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || cPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        logInfo("Direct repeat has more than "<<CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD<<" of bases the same", 6);
        return true;   
    }
    return false;
}

bool isSpacerAndDirectRepeatSimilar(DirectRepeat &dr_match)
{
    float max_length = std::max((int)dr_match.DR_Spacer.length(), dr_match.DR_Length);
    
    float edit_distance = Levensthein_distance(dr_match.DR_Sequence, dr_match.DR_Spacer);
    float similarity = 1.0 - (edit_distance/max_length);
    logInfo("similarity between spacer: "<<dr_match.DR_Spacer<<" and direct repeat: "<<dr_match.DR_Sequence<<" is: "<<similarity, 6);
    if (similarity > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        return true;
    }
    return false;
}

bool partialStarting (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq)
{
    // search in the start of the read
    logInfo("searching the beginning of the read for partials", 5);
    int tetra_start = 0;
    std::string tetramer = seq.substr(tetra_start, 4);
    // check for its presence in the DR
    size_t index = dr_match.DR_Sequence.find(tetramer);
    
    // if yes then find its maximal position ( it should reach the end of the DR)
    if (index != string::npos) 
    {
        // find max pos
        int mismatch = 0;
        // find max pos allowing for mismatches 
        while (mismatch <= CRASS_DEF_MAX_LOST_SOULS_MISMATCHES && index <= (dr_match.DR_Sequence.length() - 1) ) 
        {
            logInfo("dr pos: "<<index<<" char: "<<dr_match.DR_Sequence.at(index)<<" seq pos: "<<tetra_start<<" char: "<<seq.at(tetra_start), 9);
            
            if (dr_match.DR_Sequence.at(index) != seq.at(tetra_start)) mismatch++;
            index++;
            tetra_start++;
        }
        if (index == (dr_match.DR_Sequence.length())) 
        {
            // we have reached the end of the DR
            logInfo("reached the end of the DR!", 9);
            // tetra start will now equal the final base in the partial DR
            logInfo("spacer: "<<tetra_start<<" : "<<dr_match.DR_StartPos - 1, 9 );
            logInfo(seq.substr(tetra_start, (dr_match.DR_StartPos - tetra_start)), 9);
            
            tmp_holder->RH_StartStops.push_back(0);
            tmp_holder->RH_StartStops.push_back(tetra_start);
            
            return true;
        }
    }
    return false;
}

bool partialEnding (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq)
{
    // search in the end of the read
    logInfo("searching the end of the read for partials", 5);
    // cut a tetramer
    int tetra_start = seq.length() - 5;
    std::string tetramer = seq.substr(tetra_start, 4);
    // check for its presence in the DR
    int index = dr_match.DR_Sequence.find(tetramer);
    
    // if yes then find its maximal position ( it should reach the end of the DR)
    if (index != string::npos) 
    {
        int mismatch = 0;
        // find max pos allowing for mismatches 
        while (mismatch <= CRASS_DEF_MAX_LOST_SOULS_MISMATCHES && index >= 0 ) 
        {
            logInfo("dr pos: "<<index<<" char: "<<dr_match.DR_Sequence.at(index)<<" seq pos: "<<tetra_start<<" char: "<<seq.at(tetra_start), 10);
            if (dr_match.DR_Sequence.at(index) != seq.at(tetra_start)) mismatch++;
            index--;
            tetra_start--;
        }
        if (index < 0) 
        {
            logInfo("reached the start of the DR!", 8);
            // we have reached the start of the DR
            // tetra start will now equal the first base in the partial DR
            logInfo("spacer: "<<dr_match.DR_EndPos<<" : "<<tetra_start, 9 );
            logInfo(seq.substr(dr_match.DR_EndPos, (tetra_start - dr_match.DR_EndPos + 1)), 9);
            tmp_holder->RH_StartStops.push_back(tetra_start);
            tmp_holder->RH_StartStops.push_back(seq.length());
            return true;
        }
    }
    return false;
}

////**************************************
//// modify the positions of a crispr
////**************************************
//
//// copied from CRT source code
//
//int getActualRepeatLength(std::vector<int> &candidateCRISPR, std::string &read, int searchWindowLength, int minSpacerLength)
//{
//    int numRepeats = candidateCRISPR.size();
//    logInfo("CRT number of repeats: "<<numRepeats, 10);
//    int firstRepeatStartIndex = candidateCRISPR.front();
//    int lastRepeatStartIndex = candidateCRISPR.back();
//    
//    logInfo("first repeat start index: "<<firstRepeatStartIndex<<" last repeat start index: "<<lastRepeatStartIndex, 10);
//    
//    int shortestRepeatSpacing = candidateCRISPR.at(1) - candidateCRISPR.at(0);
//    
//    for (int i = 0; i < numRepeats - 1; i++)
//    {
//        int currRepeatIndex = candidateCRISPR.at(i);
//        int nextRepeatIndex = candidateCRISPR.at(i + 1);
//        int currRepeatSpacing = nextRepeatIndex - currRepeatIndex;
//        if (currRepeatSpacing < shortestRepeatSpacing)
//        {
//            shortestRepeatSpacing = currRepeatSpacing;
//        }
//    }
//    logInfo("shortest repeat spacing: "<<shortestRepeatSpacing, 10);
//    int sequenceLength = read.length() - 1;
//    logInfo("CRT read length: "<<sequenceLength, 10);
//    //equal to length of search string
//    int rightExtensionLength = 0;
//    
//    int currRepeatStartIndex;
//    std::string currRepeat;
//    int charCountA, charCountC, charCountT, charCountG;
//    charCountA = charCountC = charCountT = charCountG = 0;
//    float threshold;
//    bool done = false;
//    
//    threshold = .75;
//    
//    logInfo("last comparable base: "<<lastRepeatStartIndex + rightExtensionLength + searchWindowLength, 10);
//    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
//    while (!done && /*(rightExtensionLength <= maxRightExtensionLength) && */((lastRepeatStartIndex + rightExtensionLength + searchWindowLength) < sequenceLength))
//    {
//        for (int k = 0; k < candidateCRISPR.size() - 1; k++ )
//        {
//            currRepeatStartIndex = candidateCRISPR.at(k);
//            
//            char lastChar = read.at(currRepeatStartIndex + rightExtensionLength + searchWindowLength);
//            logInfo("index: "<<k<<" last char: "<<lastChar<<" at position: "<<currRepeatStartIndex + rightExtensionLength + searchWindowLength, 10);
////            currRepeat = read.substr(currRepeatStartIndex, currRepeatStartIndex + rightExtensionLength);
////            char lastChar = currRepeat.at(currRepeat.length() - 1);
//            
//            if (lastChar == 'A')   charCountA++;
//            if (lastChar == 'C')   charCountC++;
//            if (lastChar == 'T')   charCountT++;
//            if (lastChar == 'G')   charCountG++;
//        }
//        
//        double percentA = (double)charCountA/candidateCRISPR.size();
//        double percentC = (double)charCountC/candidateCRISPR.size();
//        double percentT = (double)charCountT/candidateCRISPR.size();
//        double percentG = (double)charCountG/candidateCRISPR.size();
//        
//        if ( (percentA >= threshold) || (percentC >= threshold) || (percentT >= threshold) || (percentG >= threshold) )
//        {
//            rightExtensionLength++;
//            charCountA = charCountC = charCountT = charCountG = 0;
//        }
//        else
//        {
//            done = true;
//        }
//        logInfo("right extension: "<<rightExtensionLength, 10);
//    }
//    //rightExtensionLength--;
//    
//    logInfo("right extension length: "<<rightExtensionLength, 10);
//    
//    int leftExtensionLength = 0;
//    charCountA = charCountC = charCountT = charCountG = 0;
//    done = false;
//    
//    int maxLeftExtensionLength = shortestRepeatSpacing - minSpacerLength - rightExtensionLength;
//    
//    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
//    while (!done && /*(leftExtensionLength <= maxLeftExtensionLength) && */(firstRepeatStartIndex - leftExtensionLength >= 0) )
//    {
//        for (int k = 0; k < candidateCRISPR.size() - 1; k++ )
//        {
//            currRepeatStartIndex = candidateCRISPR.at(k);
//            char firstChar = read.at(currRepeatStartIndex - leftExtensionLength);
//            logInfo("index: "<<k<<" first char: "<<firstChar<<" at position: "<<currRepeatStartIndex - leftExtensionLength, 10);
//
//            if (firstChar == 'A')    charCountA++;
//            if (firstChar == 'C')    charCountC++;
//            if (firstChar == 'T')    charCountT++;
//            if (firstChar == 'G')    charCountG++;
//        }
//        
//        double percentA = (double)charCountA/candidateCRISPR.size();
//        double percentC = (double)charCountC/candidateCRISPR.size();
//        double percentT = (double)charCountT/candidateCRISPR.size();
//        double percentG = (double)charCountG/candidateCRISPR.size();
//        
//        if ( (percentA >= threshold) || (percentC >= threshold) || (percentT >= threshold) || (percentG >= threshold) )
//        {
//            leftExtensionLength++;
//            charCountA = charCountC = charCountT = charCountG = 0;
//        }
//        else
//        {
//            done = true;
//        }
//        logInfo("left extension: "<<leftExtensionLength, 10);
//    }
//    leftExtensionLength--;
//    
//    logInfo("left extension length: "<<leftExtensionLength, 10);
//    
//    for (int m = 0; m < candidateCRISPR.size(); m++)
//    {
//        candidateCRISPR.at(m) = candidateCRISPR.at(m) - leftExtensionLength;
//    }
//    
//    return rightExtensionLength + leftExtensionLength + searchWindowLength;
//    
//}

//**************************************
// transform read to DRlowlexi
//**************************************

std::string DRLowLexi(std::string matchedRead, ReadHolder * tmp_holder)
{
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    // make sure that tere is 4 elements in the array, if not you can only cut one
    if (tmp_holder->RH_StartStops.size() == 2)
    {
        tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
        rev_comp = reverseComplement(tmp_dr);
    }
    else
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // if they both are then just take whichever is longer
        if (tmp_holder->RH_StartStops.front() == 0 && tmp_holder->RH_StartStops.back() == tmp_holder->RH_Seq.length())
        {
            int lenA = tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0);
            int lenB = tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2);
            
            if (lenA > lenB)
            {
                tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(2), (tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
        }
        // take the second
        else if (tmp_holder->RH_StartStops.front() == 0)
        {
            tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(2), (tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
        // take the first
        else 
        {
            tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
    }

    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        tmp_holder->RH_WasLowLexi = true;
        tmp_holder->RH_Seq = matchedRead;
        logInfo("DR in low lexi"<<endl<<tmp_holder->RH_Seq, 9);
        return tmp_dr;
    }
    else
    {
        tmp_holder->RH_Seq = reverseComplement(matchedRead);
        tmp_holder->reverseStartStops();
        tmp_holder->RH_WasLowLexi = false;
        logInfo("DR not in low lexi"<<endl<<tmp_holder->RH_Seq, 9);
        return rev_comp;
    }
}

void addReadHolder(ReadMap * mReads, ReadHolder * tmp_holder, std::string read_header, std::string read)
{
    logInfo("Add (header): \t" << read_header, 9);
    
    //add the header for the matched read
    tmp_holder->RH_Header = read_header;
    
    //tmp_holder->RH_StartStops
    // test drlowlexi
    std::string dr_lowlexi = DRLowLexi(read, tmp_holder);
    
    
    if (keyExists(mReads, dr_lowlexi))
    {
        // add the sequence to the map
        (*mReads)[dr_lowlexi]->push_back(tmp_holder);
    }
    else
    {
        (*mReads)[dr_lowlexi] = new ReadList();
        (*mReads)[dr_lowlexi]->push_back(tmp_holder);
    }
}

//**************************************
// lookup table shite
//**************************************

// assign a unique ID to our patterns
void addToLookup(const std::string &dr, lookupTable &patternsLookup)
{
    static int ID = 0;
    //patternsLookup[dr] = ID;
    ++ID;
}


//**************************************
// STL extensions
//**************************************

bool inline keyExists(lookupTable &patterns_hash, std::string &direct_repeat)
{
    return patterns_hash.find(direct_repeat) != patterns_hash.end();
}

bool inline keyExists(ReadMap * mReads, std::string &direct_repeat)
{
    return mReads->find(direct_repeat) != mReads->end();
}

// turn our map into a vector using just the keys
void map2Vector(lookupTable &patterns_hash, std::vector<std::string> &patterns)
{
    
    lookupTable::iterator iter = patterns_hash.begin();
    while (iter != patterns_hash.end()) 
    {
        patterns.push_back(iter->first);
        iter++;
    }
}

void map2Vector(ReadMap * mReads, std::vector<std::string> &patterns)
{

    ReadMap::reverse_iterator red_map_iter = mReads->rbegin();

    while (red_map_iter != mReads->rend()) 
    {
        patterns.push_back(red_map_iter->first);
        ++red_map_iter;
    }

}
