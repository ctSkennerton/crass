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
float bm_search_fastq_file(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &spacerLookup, lookupTable &kmerLookup, ReadMap * mReads)
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
        

        // boyer-moore search
        for (int start = 0; start < search_end; start++)
        {
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            logInfo("search begin: "<<search_begin<<" search end: "<<search_end, 8);
            
            if (search_begin >= search_end ) break;
            
            
            std::string query_word = read.substr(start, opts.lowDRsize);
            std::string subject_word = read.substr(search_begin, (opts.highDRsize + opts.highSpacerSize));
            
            logInfo("query: "<<query_word<<" subject: "<<subject_word, 10);
            
            
            int MatchStartPos = PatternMatcher::bmpSearch( subject_word, query_word );
            logInfo("bm return: "<<MatchStartPos, 8);
            
            if (MatchStartPos > -1) 
            {
                int EndPos = MatchStartPos + search_begin + opts.lowDRsize;
                int MatchEndPos = start + opts.lowDRsize;
                
                dr_match.DR_StartPos = MatchStartPos + search_begin;
                dr_match.DR_MatchStartPos = start;
                dr_match.DR_MatchEndPos = MatchEndPos;
                dr_match.DR_EndPos = EndPos;
                
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                ++EndPos;
                ++MatchEndPos;
                if (EndPos <= seq_length) 
                {
                    //read through the subsuquent bases untill they don't match
                    logInfo("read: "<<seq->name.s<<" len: "<<read.length(), 10);
                    while (read.at(MatchEndPos) == read.at(EndPos)) 
                    {
                        logInfo("match end pos: "<<MatchEndPos<<" end pos: "<<EndPos, 10);
                        logInfo(read.at(MatchEndPos) << " : " << MatchEndPos << " == " << read.at(EndPos) << " : " << EndPos, 10);

                        
                        dr_match.DR_MatchEndPos = MatchEndPos;
                        dr_match.DR_EndPos = EndPos;
                        ++EndPos;
                        ++MatchEndPos;
                        if (EndPos > seq_length) break;
                    }
                    if (!(cut_direct_repeat_sequence(dr_match, opts, read)))
                    {
                        // inctrement the for loop so that it begins at the start of the 
                        // matched kmer/direct repeat
                        // minus 1 cause it will be incremented again at the top of the for loop
                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
                        continue;
                    } 
                    else 
                    {                        
                        // create the read holder
                        ReadHolder * tmp_holder = new ReadHolder();
                        
                        addReadHolder(mReads, tmp_holder, &dr_match, read_header, read);

                        start = dr_match.DR_StartPos - 1;
                        dr_match.reset();
                        continue;
                    }
                }
            }
        }
        match_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<input_fastq, 1);
    return total_base / match_counter;
}

bool inline check_mismatch( int &temp_mismatch, const options &opts)
{
    if (temp_mismatch > opts.max_mismatches)
    {
        return false;
    }
    return true;
}

// bitap functions
// only applicable for short reads
// longer reads have enough information in them already
bool direct_repeat_init( DirectRepeat &dr_match, ReadMatch &match_info, int &temp_mismatch, const options &opts )
{
    
    dr_match.DR_NumMismatches = match_info.RM_NumMismatches;
    dr_match.DR_StartPos = match_info.RM_StartPos;
    dr_match.DR_EndPos = match_info.RM_EndPos;
    dr_match.DR_MatchStartPos = match_info.RM_MatchStartPos;
    dr_match.DR_MatchEndPos = match_info.RM_MatchEndPos;
    
    return true;  
}

bool direct_repeat_extend ( DirectRepeat &dr_match, ReadMatch &match_info, int &temp_mismatch, const options &opts)
{
    temp_mismatch += match_info.RM_NumMismatches;
    
    if (check_mismatch(temp_mismatch, opts))
    {
        dr_match.DR_NumMismatches = temp_mismatch;
        dr_match.DR_EndPos = match_info.RM_EndPos;
        dr_match.DR_MatchEndPos = match_info.RM_MatchEndPos;
        
    }
    else
    {
        return false;
    }
    
    return true;
}

bool update_word_bitap(ReadMatch &match_info, int &search_begin, int &start, const options &opts, DirectRepeat &dr, std::string &subject_word, int &temp_mismatch)
{
    match_info.RM_MatchStartPos = start;
    match_info.RM_MatchEndPos = start + opts.lowDRsize;
    match_info.RM_SubstrEnd = match_info.RM_SubstrStart + opts.lowDRsize/*strlen(opts.search_pattern)*/;
    match_info.RM_StartPos = (int) (match_info.RM_SubstrStart - subject_word.c_str()) + search_begin/*seq->seq.s)*/;
    match_info.RM_EndPos = (int) ( match_info.RM_StartPos + opts.lowDRsize/*strlen(opts.search_pattern)*/ );
    
    temp_mismatch = match_info.RM_NumMismatches;
    
    if (!(check_mismatch(temp_mismatch, opts)))
    {
        return false;
    }
    else if (dr.DR_StartPos == 0 && dr.DR_EndPos == 0)
    {
        direct_repeat_init(dr, match_info, temp_mismatch, opts);
    }
    else
    {
        if (match_info.RM_EndPos != (dr.DR_EndPos + 1)) 
        {
            return false;
        }
        direct_repeat_extend(dr, match_info, temp_mismatch, opts);
    }
    return true;
}

float bitap_search_fastq_file(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &spacerLookup, lookupTable &kmerLookup, ReadMap *mReads) 
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
        
        total_base += seq_length;
        
        dr_match.reset();
        
        std::string read = seq->seq.s;
        
        
        
        // initalize the read match with default values
        //match_init(&match_info);
        
        int temp_mismatch = 0;
        bitapType b;
        
        for (int start = 0; start < search_end; start++)
        {
            ReadMatch match_info;
            // don't search too far into the read if we don't need to
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            logInfo("search begin: "<<search_begin<<" search end: "<<search_end, 10);
            
            if (search_begin == search_end )
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
            
            if (!(NULL == (match_info.RM_SubstrEnd = FindWithBitap (&b, subject_word.c_str(), subject_word.length(), opts.max_mismatches, &match_info.RM_NumMismatches, &match_info.RM_SubstrStart)))) 
            {
                if (!(update_word_bitap(match_info, search_begin, start, opts, dr_match, subject_word, temp_mismatch)))
                {
                    if (!(cut_direct_repeat_sequence(dr_match, opts, read)))
                    {
                        dr_match.reset();
                        //++match_counter;
                        continue;
                    } 
                    else 
                    {
                        // create the read holder
                        ReadHolder * tmp_holder = new ReadHolder;
                        
                        addReadHolder(mReads, tmp_holder, &dr_match, seq->name.s, read);

                        dr_match.reset();
                    }
                }
            }
            DeleteBitap (&b);
        }
        
        // skip to the next read if the direct repeat or spacer are not the right size
        if (!(cut_direct_repeat_sequence(dr_match, opts, read)))
        {
            //++match_counter;
            continue;
        } 
        else 
        {
            // create the read holder
            ReadHolder * tmp_holder = new ReadHolder;
            
            addReadHolder(mReads, tmp_holder, &dr_match, seq->name.s, read);
        }
        ++match_counter;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    
    
    return total_base / match_counter;
}

void read_for_multimatch(const char *input_fastq, const options &opts, ReadMap *mReads )
{
    std::vector<std::string> patterns;
    
    map_to_vector(mReads, patterns);
    
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
        size_t endPos = -1;
        
        std::string read = seq->seq.s;
        
        dr_match.DR_Sequence = search.Search(strlen(seq->seq.s), seq->seq.s, patterns, endPos);
        
        dr_match.DR_StartPos = endPos;
        if (endPos != -1)
        {
            dr_match.DR_EndPos = endPos + dr_match.DR_MatchSequence.length();
            
            
            // create the read holder
            ReadHolder * tmp_holder = new ReadHolder;
            
            addReadHolder(mReads, tmp_holder, &dr_match, seq->name.s, read);
        }
    }
    logInfo("finnished multi pattern matcher", 1);
}

//**************************************
// kmer operators
//**************************************

void cutLeftKmer( string &read, int &start, int &end, lookupTable &inputLookup, const options &opts)
{    
    string left_kmer = read.substr((start - opts.kmer_size), opts.kmer_size);
    
    if (!(keyExists(inputLookup, left_kmer))) 
    {
        addToLookup(left_kmer, inputLookup);
    }
}

void cutRightKmer( string &read, int &start, int &end, lookupTable &inputLookup, const options &opts)
{
    string right_kmer = read.substr(end, opts.kmer_size);
    
    if (!(keyExists(inputLookup, right_kmer))) 
    {
        addToLookup(right_kmer, inputLookup);
    }
}

// TODO change lookuptable to proper type for spacer
void cutSpacerKmers( std::string &spacerSeq, lookupTable &spacerLookup, const options &opts)
{
    string right_kmer = spacerSeq.substr(0, opts.kmer_size);
    
    if (!(keyExists(spacerLookup, right_kmer))) 
    {
        addToLookup(right_kmer, spacerLookup);
    }
    
    string left_kmer = spacerSeq.substr((spacerSeq.length() - opts.kmer_size), opts.kmer_size);
    
    if (!(keyExists(spacerLookup, left_kmer))) 
    {
        addToLookup(left_kmer, spacerLookup);
    }
}

bool cut_direct_repeat_sequence(DirectRepeat &dr_match, const options &opts, string &read)
{
    if (!(check_dr_and_spacer_length(opts, dr_match)))
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

bool check_dr_and_spacer_length(const options &opts, DirectRepeat &dr_match)
{
    dr_match.DR_Length = dr_match.DR_EndPos - dr_match.DR_StartPos;
    int spacer_length = dr_match.DR_StartPos - dr_match.DR_MatchEndPos;
    logInfo("DR len: "<<dr_match.DR_Length<<" SP length: "<<spacer_length, 10);
    // check if the direct repeat is in the right size range
    if ( !(dr_match.DR_Length >= opts.lowDRsize && dr_match.DR_Length <= opts.highDRsize))
    {
        return false; 
    }
    
    // check if the spacer is in the right size range
    else if (!(spacer_length > opts.lowSpacerSize && spacer_length < opts.highSpacerSize))
    {
        return false;
    }
    return true; 
}

//**************************************
// transform read to DRlowlexi
//**************************************

std::string DRLowLexi(std::string matchedRead, DirectRepeat * dr_match, ReadHolder * tmp_holder)
{
    logInfo("Read Found:"<<endl<<matchedRead, 5)
    
    std::string rev_comp = reverseComplement(dr_match->DR_Sequence);
    
    if (dr_match->DR_Sequence < rev_comp/*lexicographical_compare(dr_match->DR_Sequence.begin(), dr_match->DR_Sequence.end(), rev_comp.begin(), rev_comp.end())*/)
    {
        
        // the direct repeat is in it lowest lexicographical form
        tmp_holder->RH_WasLowLexi = true;
        tmp_holder->RH_Seq = matchedRead;
        logInfo("DR in low lexi"<<endl<<tmp_holder->RH_Seq, 8);
        return dr_match->DR_Sequence;
    }
    else
    {
        tmp_holder->RH_Seq = reverseComplement(matchedRead);
        tmp_holder->RH_WasLowLexi = false;
        logInfo("DR not in low lexi"<<endl<<tmp_holder->RH_Seq, 8);
        return rev_comp;
    }

}


void addReadHolder(ReadMap * mReads, ReadHolder * tmp_holder, DirectRepeat * dr_match, std::string read_header, std::string read)
{
    logInfo(read_header, 2);
    logInfo("\tdirect repeat found:", 2);
    logInfo("\t"<<dr_match->DR_StartPos<<"\t"<<dr_match->DR_Sequence<<"\t"<<dr_match->DR_EndPos, 2);
    
    //add the header for the matched read
    tmp_holder->RH_Header = read_header;
    // test drlowlexi
    std::string dr_lowlexi = DRLowLexi(read, dr_match, tmp_holder);
    
    // populate the start stop list
    (tmp_holder->RH_StartStops).push_back( dr_match->DR_MatchStartPos);
    (tmp_holder->RH_StartStops).push_back( dr_match->DR_MatchEndPos);
    (tmp_holder->RH_StartStops).push_back(dr_match->DR_StartPos);
    (tmp_holder->RH_StartStops).push_back(dr_match->DR_EndPos);
    
    if (keyExists(mReads, dr_match->DR_Sequence))
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
// Trimming
//**************************************

// input a map of direct repeat sequence to vector of read holder objects
void drTrim(ReadMap * mReads)
{
    // load in all of the read holder objects
    
    
    // sort the read holders based on length of direct repeat
    
    
    // compare the sequence of the shorter direct repeat to the longer direct repeat
    
    
    // if it is a substring then check if the spacer contains the rest of the direct repeat
    
    
    // if yes then update the positions of the direct repeat and spacer
    
    
    // create a crispr node object and add into a node manager
    
    
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
                PRG_NAME, inputFile);
                exit(1);
    }
    
    if ( (fp == NULL) && (strcmp(inputFile, "-") == 0) ) {
        fprintf(stderr, "%s : [ERROR] Could not open stdin for reading.\n",
                PRG_NAME);
                exit(1);
    }
    return fp;
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
void map_to_vector(lookupTable &patterns_hash, std::vector<std::string> &patterns)
{
    
    lookupTable::iterator iter = patterns_hash.begin();
    while (iter != patterns_hash.end()) 
    {
        patterns.push_back(iter->first);
        iter++;
    }
}

void map_to_vector(ReadMap * mReads, std::vector<std::string> &patterns)
{

    ReadMap::reverse_iterator red_map_iter = mReads->rbegin();

    while (red_map_iter != mReads->rend()) 
    {
        patterns.push_back(red_map_iter->first);
        ++red_map_iter;
    }

}