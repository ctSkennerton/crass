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
#include "WuManber.h"
#include "libbitap.h"
#include "bm.h"
#include "kseq.h"
#include "SeqUtils.h"

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

//**************************************
// RepeatMatch
//**************************************
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
float bmSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
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
            
            logInfo("search begin: " << search_begin << " search end: " << search_end, 8);
            
            if (search_begin >= search_end ) break;
            
            std::string query_word = read.substr(start, opts.lowDRsize);
            std::string subject_word = read.substr(search_begin);
            
            logInfo("query: " << query_word << " subject: " << subject_word, 10);
            
            int MatchStartPos = PatternMatcher::bmpSearch( subject_word, query_word );
            logInfo("bm return: " << MatchStartPos, 8);
            
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

                        tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchStartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchEndPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                        match_found = true;
                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
//                        continue;
                    } 
                    else 
                    {    
                        // increment the for loop so that it begins at the start of the 
                        // matched kmer/direct repeat
                        // minus 1 cause it will be incremented again at the top of the for loop
                        match_found = false;
                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
//                        continue;
                    }
                }
            }
        }
        if (match_found)
        {
            readsFound[read_header] = true;
            addReadHolder(mReads, mStringCheck, tmp_holder, read_header, read);
        }
        else
        {
            delete tmp_holder;
        }

        match_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<input_fastq, 1);
    return total_base / match_counter;
}

float bitapSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap *mReads, StringCheck * mStringCheck) 
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

        
        int temp_mismatch = 0;
        bitapType b;
        
        ReadHolder * tmp_holder = new ReadHolder;

        for (int start = 0; start < search_end; start++)
        {
            ReadMatch match_info;
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
            
            if (NULL != (match_info.RM_SubstrEnd = FindWithBitap (&b, subject_word.c_str(), subject_word.length(), opts.max_mismatches, &match_info.RM_NumMismatches, &match_info.RM_SubstrStart))) 
            {
                dr_match.DR_MatchStartPos = start;
                dr_match.DR_StartPos = (int) (match_info.RM_SubstrStart - subject_word.c_str()) + search_begin;

                
                match_info.RM_MatchEndPos = start + opts.lowDRsize;
                match_info.RM_SubstrEnd = match_info.RM_SubstrStart + opts.lowDRsize;
                match_info.RM_EndPos = (int) ( dr_match.DR_StartPos + opts.lowDRsize );
                
                            
                ++match_info.RM_MatchEndPos;
                ++match_info.RM_EndPos;
                
                while (match_info.RM_NumMismatches <= opts.max_mismatches && match_info.RM_EndPos < seq_length) 
                {
                    char dr_A = read.at(match_info.RM_MatchEndPos);
                    char dr_B = read.at(match_info.RM_EndPos);
                    
                    logInfo("Match end pos: "<<match_info.RM_MatchEndPos<<" end pos: "<<match_info.RM_EndPos, 10);
                    logInfo(read.at(match_info.RM_MatchEndPos) << " : " << match_info.RM_MatchEndPos << " == " << read.at(match_info.RM_EndPos) << " : " << match_info.RM_EndPos, 10);
                    
                    if (dr_A != dr_B)
                    {
                        match_info.RM_NumMismatches++;
                        if (match_info.RM_NumMismatches <= opts.max_mismatches) break;
                    }
                    
                    ++match_info.RM_MatchEndPos;
                    ++match_info.RM_EndPos;
                }
                
                dr_match.DR_MatchEndPos = match_info.RM_MatchEndPos;
                dr_match.DR_NumMismatches = match_info.RM_NumMismatches;
                dr_match.DR_EndPos = match_info.RM_EndPos;
                
                if (cutDirectRepeatSequence(dr_match, opts, read))
                {

                    patterns_hash[dr_match.DR_MatchSequence] = true;
                    
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchStartPos);
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchEndPos);
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
                    // TODO True or false?
                    match_found = false;
                    dr_match.reset();
                    continue;
                }
            }
            DeleteBitap (&b);
        }
        
        if (match_found)
        {
            readsFound[seq->name.s] = true;
            addReadHolder(mReads, mStringCheck, tmp_holder, seq->name.s, read);
        }
        else
        {
            delete tmp_holder;
        }
        match_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    
    std::cout<<total_base<<" "<<match_counter<<endl;
    
    return total_base / match_counter;
}

void scanForMultiMatches(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    std::vector<std::string> patterns;
    
    map2Vector(patterns_hash, patterns);
    
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
    
    int l;
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
    
        DirectRepeat dr_match;        
        
        //initialize with an impossible number
        int start_pos = -1;
        
        std::string read = seq->seq.s;

        dr_match.DR_Sequence = search.Search(strlen(seq->seq.s), seq->seq.s, patterns, start_pos);
        dr_match.DR_StartPos = start_pos;
        
        if (start_pos != -1)
        {
            std::string header = seq->name.s;
            if (readsFound.find(header) == readsFound.end())
            //if (!(keyExists(readsFound, header)))
            {
                dr_match.DR_EndPos = start_pos + dr_match.DR_Sequence.length();
                
                // TODO: change this to some smarter logic
                // really silly way or breaking it up!!
                int thirds = read.length()/3;
                int f_start = -1;
                int f_end = -1;
                
                // create the read holder
                ReadHolder * tmp_holder = new ReadHolder;
                
                if (dr_match.DR_EndPos <= thirds)
                {
                    // first third
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                    tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                    if(partialEnding(dr_match, tmp_holder, read, &f_start, &f_end))
                    {
                        int dist = f_start - dr_match.DR_EndPos;
                        if ((dist >= opts.lowSpacerSize) and (dist <= opts.highSpacerSize))
                        {
                            tmp_holder->RH_StartStops.push_back(f_start);
                            tmp_holder->RH_StartStops.push_back(f_end);
                        }
                    }
                }
                else if (dr_match.DR_EndPos >= thirds*2)
                {
                    // last third
                    if(partialStarting(dr_match, tmp_holder, read, &f_start, &f_end))
                    {
                        int dist = dr_match.DR_StartPos - f_end;
                        if ((dist >= opts.lowSpacerSize) and (dist <= opts.highSpacerSize))
                        {
                            tmp_holder->RH_StartStops.push_back(f_start);
                            tmp_holder->RH_StartStops.push_back(f_end);
                        }
                    }
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_StartPos);
                    tmp_holder->RH_StartStops.push_back( dr_match.DR_EndPos);
                }
                else
                {
                    // middle third
                    if (partialStarting(dr_match, tmp_holder, read, &f_start, &f_end))
                    {
                        int dist = dr_match.DR_StartPos - f_end;
                        if ((dist >= opts.lowSpacerSize) and (dist <= opts.highSpacerSize))
                        {
                            tmp_holder->RH_StartStops.push_back(f_start);
                            tmp_holder->RH_StartStops.push_back(f_end);
                        }
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                    }
                    else
                    {
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                        if(partialEnding(dr_match, tmp_holder, read, &f_start, &f_end))
                        {
                            int dist = f_start - dr_match.DR_EndPos;
                            if ((dist >= opts.lowSpacerSize) and (dist <= opts.highSpacerSize))
                            {
                                tmp_holder->RH_StartStops.push_back(f_start);
                                tmp_holder->RH_StartStops.push_back(f_end);
                            }
                        }
                    }
                }
                addReadHolder(mReads, mStringCheck, tmp_holder, seq->name.s, read);
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
    
    if (!(checkDRAndSpacerLength(opts, dr_match)) || isLowComplexity(dr_match, read))
    {
        return false;
    }

    // if the length of both spacer and direct repeat are okay cut the subsequences
    else
    {
        dr_match.DR_Sequence = read.substr(dr_match.DR_StartPos, (dr_match.DR_EndPos - dr_match.DR_StartPos + 1));
        dr_match.DR_MatchSequence = read.substr(dr_match.DR_MatchStartPos, (dr_match.DR_MatchEndPos - dr_match.DR_MatchStartPos + 1));
        dr_match.DR_Spacer = read.substr(dr_match.DR_MatchEndPos, (dr_match.DR_StartPos - dr_match.DR_MatchEndPos));
    } 
    return true;
}

bool checkDRAndSpacerLength(const options &opts, DirectRepeat &dr_match)
{
    logInfo("dr end pos: "<<dr_match.DR_EndPos<<" dr start pos: "<<dr_match.DR_StartPos, 10);
    
    int spacer_length = dr_match.DR_StartPos - dr_match.DR_MatchEndPos;
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

bool isLowComplexity(DirectRepeat &dr_match, std::string read)
{
    int cCount = 0;
    int gCount = 0;
    int aCount = 0;
    int tCount = 0;
    
    float aPercent;
    float cPercent;
    float gPercetn;
    float tPercent;

    
    int i = dr_match.DR_StartPos;
    while (i <= dr_match.DR_EndPos) 
    {
        switch (read[i]) 
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
        i++;
    }
    aPercent = aCount/dr_match.DR_Length;
    tPercent = tCount/dr_match.DR_Length;
    gPercetn = gCount/dr_match.DR_Length;
    cPercent = cCount/dr_match.DR_Length;
    
    if (aPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || tPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || gPercetn > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || cPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        return true;   
    }
    return false;
}

bool partialStarting (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq, int * f_start, int * f_end)
{
return false;
    // search in the start of the read
    logInfo("searching the beginning of the read for partials", 5);
    int tetra_start = 0;
    std::string tetramer = seq.substr(tetra_start, 4);
    // check for its presence in the DR
    int index = dr_match.DR_Sequence.find(tetramer);
    
    // if yes then find its maximal position ( it should reach the end of the DR)
    if (index != string::npos) 
    {
        // find max pos
        int mismatch = 0;
        // find max pos allowing for mismatches 
        while (mismatch <= CRASS_DEF_MAX_LOST_SOULS_MISMATCHES && index < (dr_match.DR_Sequence).length()) 
        {
            logInfo("dr pos: "<<index<<" char: "<<dr_match.DR_Sequence.at(index)<<" seq pos: "<<tetra_start<<" char: "<<seq.at(tetra_start), 9);
            
            if (dr_match.DR_Sequence.at(index) != seq.at(tetra_start)) mismatch++;
            index++;
            tetra_start++;
        }
        if ((index == ((dr_match.DR_Sequence).length())) && (mismatch <= CRASS_DEF_MAX_LOST_SOULS_MISMATCHES)) 
        {
            // we have reached the end of the DR
            logInfo("reached the end of the DR!", 9);
            // tetra start will now equal the final base in the partial DR
            logInfo("spacer: "<<tetra_start<<" : "<<dr_match.DR_StartPos - 1, 9 );
            logInfo(seq.substr(tetra_start, (dr_match.DR_StartPos - tetra_start)), 9);
            
            *f_start = 0;
            *f_end = tetra_start;
            return true;
        }
    }
    return false;
}

bool partialEnding (DirectRepeat &dr_match, ReadHolder *tmp_holder, std::string &seq, int * f_start, int * f_end)
{
return false;
    // search in the end of the read
    logInfo("searching the end of the read for partials", 5);
    // cut a tetramer
    int tetra_start = seq.length() - 5;
    std::string tetramer = seq.substr(tetra_start, 4);
    // check for its presence in the DR
    int index = dr_match.DR_Sequence.rfind(tetramer);

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
        if (index == -1 && (mismatch <= CRASS_DEF_MAX_LOST_SOULS_MISMATCHES)) 
        {
            logInfo("reached the start of the DR!", 8);
            // we have reached the start of the DR
            // tetra start will now equal the first base in the partial DR
            logInfo("spacer: "<<dr_match.DR_EndPos<<" : "<<tetra_start, 9 );
            logInfo(seq.substr(dr_match.DR_EndPos, (tetra_start - dr_match.DR_EndPos + 1)), 9);
            
            *f_start = tetra_start;
            *f_end = seq.length() - 1;
            return true;
        }
    }
    return false;
}

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

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmp_holder, std::string read_header, std::string read)
{
    logInfo("Add (header): \t" << read_header, 9);
    
    //add the header for the matched read
    tmp_holder->RH_Header = read_header;
    std::string dr_lowlexi = DRLowLexi(read, tmp_holder);
    StringToken st = mStringCheck->getToken(dr_lowlexi);
    if(0 == st)
    {
        // new guy
        st = mStringCheck->addString(dr_lowlexi);
        (*mReads)[st] = new ReadList();
    }
    (*mReads)[st]->push_back(tmp_holder);
}

//**************************************
// system
//**************************************

gzFile getFileHandle(const char * inputFile)
{
    gzFile fp;
    if ( strcmp(inputFile, "-") == 0 ) {
        fp = gzdopen(fileno(stdin), "r");
    }
    else {
        fp = gzopen(inputFile, "r");
    }
    
    if ( (fp == NULL) && (strcmp(inputFile, "-") != 0) ) {
        fprintf(stderr, "%s : [ERROR] Could not open FASTQ '%s' for reading.\n",
                CRASS_DEF_PRG_NAME, inputFile);
                exit(1);
    }
    
    if ( (fp == NULL) && (strcmp(inputFile, "-") == 0) ) {
        fprintf(stderr, "%s : [ERROR] Could not open stdin for reading.\n",
                CRASS_DEF_PRG_NAME);
                exit(1);
    }
    return fp;
}

//**************************************
// STL extensions
//**************************************

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

