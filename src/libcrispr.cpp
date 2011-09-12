// File: libcrispr.cpp
// Original Author: Michael Imelfort 2011  :)
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
#include <sstream>
#include <zlib.h>  
#include <fstream>
#include <cmath>
#include <fcntl.h>
#include <stdlib.h>


// local includes
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "WuManber.h"
#include "PatternMatcher.h"
#include "SeqUtils.h"
#include "Levensthein.h"
#include "kseq.h"
#include "StlExt.h"


/* 
declare the type of file handler and the read() function
as described here:
http://lh3lh3.users.sourceforge.net/parsefastq.shtml

THIS JUST DEFINES A BUNCH OF **templated** structs

*/
    KSEQ_INIT(gzFile, gzread)

#pragma mark Search Functions
READ_TYPE decideWhichSearch(const char *inputFastq, float * aveReadLength)
{
    //-----
    // Wrapper used for searching reads for DRs
    // depending on the length of the read. this funciton may use the boyer moore algorithm
    // of the CRT search algorithm
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;

    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        int seq_length = (int)strlen(seq->seq.s);
        total_base += seq_length;
        read_counter++;
        if(read_counter > 100)
            break;
    }
    gzclose(fp);
    
    *aveReadLength = total_base / read_counter;
    logInfo("Average read length (of the first 100 reads): "<< *aveReadLength, 2);
    if((total_base / read_counter) > CRASS_DEF_READ_LENGTH_CUTOFF)
    {
        return LONG_READ;
    }
    return SHORT_READ;
}


// CRT search

void longReadSearch(const char *inputFastq, const options& opts, ReadMap * mReads, StringCheck * mStringCheck)
{
    //-----
    // Code lifted from CRT, ported by connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        // grab a readholder
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s,seq->name.s);
        bool match_found = false;

        if (opts.removeHomopolymers)
        {
                // RLE is necessary...
                tmp_holder->encode();
        }
        std::string read = tmp_holder->getSeq();
        
        // get the length of this sequence
        unsigned int seq_length = (unsigned int)read.length();
        
        total_base += seq_length;

        //the mumber of bases that can be skipped while we still guarantee that the entire search
        //window will at some point in its iteration thru the sequence will not miss a any repeat
        unsigned int skips = opts.lowDRsize - (2 * opts.searchWindowLength - 1);
        if (skips < 1)
        {
            skips = 1;
        }

        int searchEnd = seq_length - opts.highDRsize - opts.highSpacerSize - opts.searchWindowLength - 1;
        
        if (searchEnd < 0) 
        {
            logWarn("Read: "<<tmp_holder->getHeader()<<" length is less than "<<opts.highDRsize + opts.highSpacerSize + opts.searchWindowLength + 1<<"bp",4);
//-DDEBUG#ifdef DEBUG 
            logWarn("This may be due to the input file containing a proportion of reads less than "<<CRASS_DEF_READ_LENGTH_CUTOFF, 4);
            logWarn(PACKAGE_NAME<<" chooses a search algorithm baised on an average length of reads",4);
            logWarn("in the beginning of the input file, not on a per read basis.", 4);
            logWarn("This read will be skipped!", 4);
//-DDEBUG#endif 
            continue;
        }
        
        for (unsigned int j = 0; j <= (unsigned int)searchEnd; j = j + skips)
        {
                        
            unsigned int beginSearch = j + opts.lowDRsize + opts.lowSpacerSize;
            unsigned int endSearch = j + opts.highDRsize + opts.highSpacerSize + opts.searchWindowLength;
            
            if (endSearch >= seq_length)
            {
                endSearch = seq_length - 1;
            }
            //should never occur
            if (endSearch < beginSearch)
            {
                endSearch = beginSearch;
            }
            
            std::string text = read.substr(beginSearch, (endSearch - beginSearch));
            std::string pattern = read.substr(j, opts.searchWindowLength);
            
            //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
            int pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);
            if (pattern_in_text_index >= 0)
            {
                tmp_holder->startStopsAdd(j,  j + opts.searchWindowLength);
                unsigned int found_pattern_start_index = beginSearch + (unsigned int)pattern_in_text_index;
                
                tmp_holder->startStopsAdd(found_pattern_start_index, found_pattern_start_index + opts.searchWindowLength);
                
                scanRight(tmp_holder, pattern, opts.lowSpacerSize, 24);
            }

            if ( (tmp_holder->numRepeats() > opts.minNumRepeats) ) //numRepeats is half the size of the StartStopList
            {
//-DDEBUG#ifdef DEBUG
                logInfo(tmp_holder->getHeader(), 8);
                logInfo("\tPassed test 1. More than "<<opts.minNumRepeats<< " ("<<tmp_holder->numRepeats()<<") repeated kmers found", 8);
//-DDEBUG#endif
                unsigned int actual_repeat_length = extendPreRepeat(tmp_holder, opts.searchWindowLength, opts.lowSpacerSize);

                if ( (actual_repeat_length >= opts.lowDRsize) && (actual_repeat_length <= opts.highDRsize) )
                {
//-DDEBUG#ifdef DEBUG

                    logInfo("\tPassed test 2. The repeat length is between "<<opts.lowDRsize<<" and "<<opts.highDRsize, 8);
//-DDEBUG#endif

                    // drop partials
                    tmp_holder->dropPartials();
                    if (qcFoundRepeats(tmp_holder))
                    {
//-DDEBUG#ifdef DEBUG
                        logInfo("Passed all tests!", 8);
                        logInfo("Potential CRISPR containing read found: "<<tmp_holder->getHeader(), 7);
                        logInfo(tmp_holder->getSeq(), 9);
                        logInfo("-------------------", 7)
//-DDEBUG#endif                            
                            match_found = true;
                            addReadHolder(mReads, mStringCheck, tmp_holder);
                            break;
                    }
                }
//-DDEBUG#ifdef DEBUG                
                else
                {
                    logInfo("\tFailed test 2. Repeat length: "<<tmp_holder->getRepeatLength() << " : " << match_found, 8); 
                }
//-DDEBUG#endif
                j = tmp_holder->back() - 1;
            }
            tmp_holder->clearStartStops();
        }
        if (!match_found) 
        {
            delete tmp_holder;
        }
        read_counter++;
    }
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    
    logInfo("finished processing file:"<<inputFastq, 1);    
    logInfo("So far " << mReads->size()<<" direct repeat variants have been found from " << read_counter << " reads", 2);

}

// boyer moore functions
void shortReadSearch(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
        
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        bool match_found = false;

        // create the read holder
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
            
        if (opts.removeHomopolymers)
        {
            // RLE is necessary...
            tmp_holder->encode();
        }
        
        std::string read = tmp_holder->getSeq();
        

        unsigned int seq_length = (unsigned int)read.length() - 1;
        unsigned int search_end = seq_length - opts.lowDRsize;
        
        total_base += seq_length;
        
        
        // boyer-moore search
        for (unsigned int first_start = 0; first_start < search_end; first_start++)
        {
            unsigned int search_begin = first_start + opts.lowDRsize + opts.lowSpacerSize;
            
            if (search_begin >= search_end ) break;
            
            int second_start = PatternMatcher::bmpSearch( read.substr(search_begin), read.substr(first_start, opts.lowDRsize) );
            
            if (second_start > -1) 
            {
                second_start += search_begin;
                unsigned int second_end = (unsigned int)second_start + opts.lowDRsize;
                unsigned int first_end = first_start + opts.lowDRsize;
                

                unsigned int next_index = second_end + 1;
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                if ( next_index <= seq_length) 
                {
                    // read through the subsuquent bases untill they don't match
                    unsigned int extenstion_length = 0;
                    while (read.at(first_end + extenstion_length) == read.at(second_end + extenstion_length)) 
                    {
                        logInfo(read.at(first_end + extenstion_length)<<" == "<<read.at(second_end + extenstion_length),8);
                        extenstion_length++;
                        next_index++;
                        if (next_index > seq_length) break;

                    }
                    logInfo(second_end << " : "<<second_start<<" : "<<extenstion_length, 8);
                    tmp_holder->startStopsAdd(first_start, first_end + extenstion_length);
                    tmp_holder->startStopsAdd(second_start, second_end + extenstion_length);
                    tmp_holder->setRepeatLength(second_end - second_start + extenstion_length);
                }
                else
                {
                    tmp_holder->startStopsAdd(first_start, first_end);
                    tmp_holder->startStopsAdd(second_start, second_end);
                    tmp_holder->setRepeatLength(second_end - second_start);
                }
                
                // the low side will always be true since we search for the lowDRSize
                if ( tmp_holder->getRepeatLength() <= opts.highDRsize )
                {

                    if ((tmp_holder->getAverageSpacerLength() >= opts.lowSpacerSize) && (tmp_holder->getAverageSpacerLength() <= opts.highSpacerSize))
                    {
                        if (qcFoundRepeats(tmp_holder)) 
                        {

                            match_found = true;
//-DDEBUG#ifdef DEBUG
                            logInfo("Potential CRISPR containing read found: "<<tmp_holder->getHeader(), 7);
                            logInfo(read, 9);
                            logInfo("-------------------", 7)
//-DDEBUG#endif
                            patternsHash[tmp_holder->repeatStringAt(0)] = true;
                            readsFound[tmp_holder->getHeader()] = true;
                            addReadHolder(mReads, mStringCheck, tmp_holder);
                            break;
                        }
                    }
//-DDEBUG#ifdef DEBUG
                    else
                    {
                        logInfo("\tFailed test 2. The spacer length is not between "<<opts.lowSpacerSize<<" and "<<opts.highSpacerSize<<": "<<tmp_holder->getAverageSpacerLength(), 8);
                    }
//-DDBUG#endif
                }
//-DDEBUG#ifdef DEBUG
                else
                {
                    logInfo("\tFailed test 1. The repeat length is larger than "<<opts.highDRsize<<": " << tmp_holder->getRepeatLength(), 8);
                }
//-DDBUG#endif
                first_start = tmp_holder->back();
            }
        }
        if (!match_found)
        {
            delete tmp_holder;
        }
        read_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("Finished processing file:"<<inputFastq, 1);
    logInfo("So far " << mReads->size()<<" direct repeat variants have been found", 2);
}


void findSingletons(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    std::vector<std::string> patterns;
    int old_number = (int)mReads->size();
    mapToVector(patternsHash, patterns);
    try {
        if (patterns.empty())
        {
            logError("No patterns in vector for multimatch");
            throw "No patterns in vector for multimatch";
        }
    } catch (char * c) {
        std::cerr<<c<<std::endl;
        return;
    }   

    
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    
    seq = kseq_init(fp);

    WuManber search;
    search.Initialize(patterns);
    
    int l;
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
        
        if (opts.removeHomopolymers) 
        {
            tmp_holder->encode();
        }
        std::string read = tmp_holder->getSeq();
        
        //initialize with an impossible number
        int start_pos = -1;
        std::string found_repeat = search.Search(read.length(), read.c_str(), patterns, start_pos);
        
        
        if (start_pos != -1)
        {
            if (readsFound.find(tmp_holder->getHeader()) == readsFound.end())
            {
//-DDEBUG#ifdef DEBUG
                logInfo("new read recruited: "<<tmp_holder->getHeader(), 7);
                logInfo(tmp_holder->getSeq(), 8);
//-DDEBUG#endif
                unsigned int DR_end = (unsigned int)start_pos + (unsigned int)found_repeat.length();
                if(DR_end >= (unsigned int)read.length())
                {
                    DR_end = (unsigned int)read.length() - 1;
                }
                tmp_holder->startStopsAdd(start_pos, DR_end);
                addReadHolder(mReads, mStringCheck, tmp_holder);
            }
        }
        else
        {
            delete tmp_holder;
        }
    }
    logInfo("Finished second iteration. An extra " << mReads->size() - old_number<<" variants were recruited", 2);
}

#pragma mark -
#pragma mark Repeat QC

int scanRight(ReadHolder * tmp_holder, std::string& pattern, unsigned int minSpacerLength, unsigned int scanRange)
{
    //-DDEBUG#ifdef DEBUG
    logInfo("Scanning Right for more repeats:", 9);
    //-DDEBUG#endif
    unsigned int start_stops_size = tmp_holder->getStartStopListSize();
    
    unsigned int pattern_length = (unsigned int)pattern.length();
    
    // final start index
    unsigned int last_repeat_index = tmp_holder->getRepeatAt(start_stops_size - 2);
    
    //second to final start index
    unsigned int second_last_repeat_index = tmp_holder->getRepeatAt(start_stops_size - 4);
    
    unsigned int repeat_spacing = last_repeat_index - second_last_repeat_index;
    
    //-DDEBUG#ifdef DEBUG
    logInfo(start_stops_size<<" : "<<pattern_length<<" : "<<last_repeat_index<<" : "<<second_last_repeat_index<<" : "<<repeat_spacing, 9);
    //-DDEBUG#endif
    
    int candidate_repeat_index, position;
    
    unsigned int begin_search, end_search;
    
    unsigned int read_length = (unsigned int)tmp_holder->getSeqLength();
    bool more_to_search = true;
    while (more_to_search)
    {
        candidate_repeat_index = last_repeat_index + repeat_spacing;
        begin_search = candidate_repeat_index - scanRange;
        end_search = candidate_repeat_index + pattern_length + scanRange;
        //-DDEBUG#ifdef DEBUG
        logInfo(candidate_repeat_index<<" : "<<begin_search<<" : "<<end_search, 9);
        //-DDEBUG#endif
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right
        unsigned int scanRightMinBegin = last_repeat_index + pattern_length + minSpacerLength;
        
        if (begin_search < scanRightMinBegin)
        {
            begin_search = scanRightMinBegin;
        }
        if (begin_search > read_length - 1)
        {
            //-DDEBUG#ifdef DEBUG
            logInfo("returning... "<<begin_search<<" > "<<read_length - 1, 9);
            //-DDEBUG#endif
            return read_length - 1;
        }
        if (end_search > read_length)
        {
            end_search = read_length;
        }
        
        if ( begin_search >= end_search)
        {
            //-DDEBUG#ifdef DEBUG
            logInfo("Returning... "<<begin_search<<" >= "<<end_search, 9);
            //-DDEBUG#endif
            return end_search;
        }
        /******************** end range checks ********************/
        
        std::string text = tmp_holder->substr(begin_search, (end_search - begin_search));
        
        //-DDEBUG#ifdef DEBUG
        logInfo(pattern<<" : "<<text, 9);
        //-DDEBUG#endif
        position = PatternMatcher::bmpSearch(text, pattern);
        
        
        if (position >= 0)
        {
            tmp_holder->startStopsAdd(begin_search + position, begin_search + position + pattern_length);
            second_last_repeat_index = last_repeat_index;
            last_repeat_index = begin_search + position;
            repeat_spacing = last_repeat_index - second_last_repeat_index;
            if (repeat_spacing < (minSpacerLength + pattern_length))
            {
                more_to_search = false;
            }
        }
        else
        {
            more_to_search = false;
        }
    }
    
    return begin_search + position;
}


unsigned int extendPreRepeat(ReadHolder * tmp_holder, int searchWindowLength, int minSpacerLength)
{
//-DDEBUG#ifdef DEBUG
    logInfo("Extending Prerepeat...", 9);
//-DDEBUG#endif
    //-----
    // Extend a preliminary repeat - return the final repeat size
    //
    // LOOK AT ME
    //!!!!!!!!!!!!!!!
    // the number of repeats
    // equal to half the size of the start stop list
    //!!!!!!!!!!!!!!
    unsigned int num_repeats = tmp_holder->numRepeats();
    tmp_holder->setRepeatLength(searchWindowLength);
    int cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * num_repeats);
    
    // make sure that we don't go below 2
    if (2 > cut_off) 
    {
        cut_off = 2;
    }
//-DDEBUG#ifdef DEBUG
    logInfo("cutoff: "<<cut_off, 9);
//-DDEBUG#endif
    
    
    // the index in the read of the first DR kmer
    unsigned int first_repeat_start_index = tmp_holder->getFirstRepeatStart();
    
    // the index in the read of the last DR kmer
    unsigned int last_repeat_start_index = tmp_holder->getLastRepeatStart();
    
    // the length between the first two DR kmers
    unsigned int shortest_repeat_spacing = tmp_holder->startStopsAt(2) - tmp_holder->startStopsAt(0);
    // loop througth all remaining members of mRepeats
    unsigned int end_index = tmp_holder->getStartStopListSize();
    
    for (unsigned int i = 4; i < end_index; i+=2)
    {
        
        // get the repeat spacing of this pair of DR kmers
        unsigned int curr_repeat_spacing = tmp_holder->startStopsAt(i) - tmp_holder->startStopsAt(i - 2);
//-DDEBUG#ifdef DEBUG
        logInfo(i<<" : "<<curr_repeat_spacing, 10);
//-DDEBUG#endif
        
        // if it is shorter than what we already have, make it the shortest
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }
//-DDEBUG#ifdef DEBUG
    logInfo("shortest repeat spacing: "<<shortest_repeat_spacing, 9);
//-DDEBUG#endif
    unsigned int right_extension_length = 0;
    // don't search too far  
    unsigned int max_right_extension_length = shortest_repeat_spacing - minSpacerLength;
//-DDEBUG#ifdef DEBUG
    logInfo("max ringt extension length: "<<max_right_extension_length, 9);
//-DDEBUG#endif
    // Sometimes we shouldn't use the far right DR. (it may lie too close to the end)
    unsigned int DR_index_end = end_index;
    unsigned int dist_to_end = tmp_holder->getSeqLength() - last_repeat_start_index - 1;
    if(dist_to_end < max_right_extension_length)
    {
//-DDEBUG#ifdef DEBUG
        logInfo("removing end partial: "<<dist_to_end<<" < "<<max_right_extension_length, 9);
//-DDEBUG#endif
        DR_index_end -= 2;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off) 
        {
            cut_off =2;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("new cutoff: "<<cut_off, 9);
//-DDEBUG#endif
    }
    std::string curr_repeat;
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
    while (max_right_extension_length > 0)
    {
        for (unsigned int k = 0; k < DR_index_end; k+=2 )
        {
//-DDEBUG#ifdef DEBUG
            logInfo(k<<" : "<<tmp_holder->getRepeatAt(k) + tmp_holder->getRepeatLength(), 10);
//-DDEBUG#endif
            // look at the character just past the end of the last repeat
            switch( tmp_holder->getSeqCharAt(tmp_holder->getRepeatAt(k) + tmp_holder->getRepeatLength()))
            {
                case 'A':
                    char_count_A++;
                    break;
                case 'C':
                    char_count_C++;
                    break;
                case 'G':
                    char_count_G++;
                    break;
                case 'T':
                    char_count_T++;
                    break;
            }
        }
        
//-DDEBUG#ifdef DEBUG
        logInfo("R: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << tmp_holder->getRepeatLength() << " : " << max_right_extension_length, 9);
//-DDEBUG#endif
        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            tmp_holder->incrementRepeatLength();
            max_right_extension_length--;
            right_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }
    
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    
    // again, not too far
    unsigned int left_extension_length = 0;
    unsigned int max_left_extension_length = shortest_repeat_spacing - minSpacerLength - tmp_holder->getRepeatLength();
    
//-DDEBUG#ifdef DEBUG
    logInfo("max left extension length: "<<max_left_extension_length, 9);
//-DDEBUG#endif
    // and again, we may not wat to use the first DR
    unsigned int DR_index_start = 0;
    if(max_left_extension_length > first_repeat_start_index)
    {
//-DDEBUG#ifdef DEBUG
        logInfo("removing start partial: "<<max_left_extension_length<<" > "<<first_repeat_start_index, 9);
//-DDEBUG#endif
        DR_index_start+=2;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off) 
        {
            cut_off = 2;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("new cutoff: "<<cut_off, 9);
//-DDEBUG#endif
    }
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (left_extension_length < max_left_extension_length)
    {
        for (unsigned int k = DR_index_start; k < end_index; k+=2 )
        {
//-DDEBUG#ifdef DEBUG
            logInfo(k<<" : "<<tmp_holder->getRepeatAt(k) - left_extension_length - 1, 10);
//-DDEBUG#endif
            switch(tmp_holder->getSeqCharAt(tmp_holder->getRepeatAt(k) - left_extension_length - 1))
            {
                case 'A':
                    char_count_A++;
                    break;
                case 'C':
                    char_count_C++;
                    break;
                case 'G':
                    char_count_G++;
                    break;
                case 'T':
                    char_count_T++;
                    break;
            }
        }
//-DDEBUG#ifdef DEBUG
        logInfo("L:" << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << tmp_holder->getRepeatLength() << " : " << left_extension_length, 9);
//-DDEBUG#endif
        
        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            tmp_holder->incrementRepeatLength();
            left_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }
    StartStopListIterator repeat_iter = tmp_holder->begin();
//-DDEBUG#ifdef DEBUG    
    logInfo("Repeat positions:", 9);
//-DDEBUG#endif
    while (repeat_iter < tmp_holder->end()) 
    {
        if(*repeat_iter < (unsigned int)left_extension_length)
        {
            *repeat_iter = 0;
            *(repeat_iter + 1) += right_extension_length;
        }
        else
        {
            *repeat_iter -= left_extension_length;
            *(repeat_iter+1) += right_extension_length;
        }
//-DDEBUG#ifdef DEBUG    
        logInfo(*repeat_iter<<","<<*(repeat_iter+1), 9);
//-DDEBUG#endif
        repeat_iter += 2;
    }
    
    return (unsigned int)tmp_holder->getRepeatLength();
    
}


//need at least two elements
bool qcFoundRepeats(ReadHolder * tmp_holder)
{
    try {
        if (tmp_holder->numRepeats() < 2) 
        {
            logError("The vector holding the repeat indexes has less than 2 repeats!");
            throw "The vector holding the repeat indexes has less than 2 repeats!";
        }
    } catch (char * c) {
        std::cerr<<c<<std::endl;
    }
    
    std::string repeat = tmp_holder->repeatStringAt(0);
    
    if (isRepeatLowComplexity(repeat)) 
    {
//-DDEBUG#ifdef DEBUG
        logInfo("\tFailed test 3. The repeat is low complexity", 8);
//-DDEBUG#endif
        return false;
    }
    
//-DDEBUG#ifdef DEBUG
    logInfo("\tPassed test 3. The repeat is not low complexity", 8);
//-DDEBUG#endif
    // test for a long or short read
    if (2 <= tmp_holder->numSpacers()) 
    {
        std::vector<std::string> spacer_vec = tmp_holder->getAllSpacerStrings();
        
        int spacer_vec_size = (int)spacer_vec.size();
        int sum_spacer_to_spacer_len_diff = 0;
        int sum_repeat_to_spacer_len_diff = 0;
        
        float sum_spacer_to_spacer_difference = 0.0;
        float sum_repeat_to_spacer_difference = 0.0;
        
        std::vector<std::string>::iterator spacer_iter = spacer_vec.begin();
        while (spacer_iter != spacer_vec.end() - 1) 
        {
            
            sum_repeat_to_spacer_difference += PatternMatcher::getStringSimilarity(repeat, *spacer_iter);
            sum_spacer_to_spacer_difference += PatternMatcher::getStringSimilarity(*spacer_iter, *(spacer_iter + 1));
            
            sum_spacer_to_spacer_len_diff += spacer_iter->size() - (spacer_iter + 1)->size();
            sum_repeat_to_spacer_len_diff += repeat.size() - spacer_iter->size();
            
            spacer_iter++;
        }
        float diff_cutoff = (spacer_vec_size * CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY);
        
        //float spacer_vec_size_as_float = (float)spacer_vec.size();
        if (sum_spacer_to_spacer_difference > diff_cutoff) 
        {
//-DDEBUG#ifdef DEBUG
            logInfo("\tFailed test 4a. Spacers are too similar: "<<sum_spacer_to_spacer_difference/spacer_vec_size<<" > "<<diff_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("\tPassed test 4a. Spacers are not too similar: "<<sum_spacer_to_spacer_difference/spacer_vec_size<<" < "<<diff_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif    
        if (sum_repeat_to_spacer_difference > diff_cutoff)
        {
//-DDEBUG#ifdef DEBUG
            logInfo("\tFailed test 4b. Spacers are too similar to the repeat: "<<sum_repeat_to_spacer_difference/spacer_vec_size<<" > "<<diff_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("\tPassed test 4b. Spacers are not too similar to the repeat: "<<sum_repeat_to_spacer_difference/spacer_vec_size<<" < "<<diff_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif    
        
        int spacer_len_cutoff = spacer_vec_size * CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF;
        if (abs((int)sum_spacer_to_spacer_len_diff) > spacer_len_cutoff) 
        {
//-DDEBUG#ifdef DEBUG 
            logInfo("\tFailed test 5a. Spacer lengths differ too much: "<<abs(sum_spacer_to_spacer_len_diff)/spacer_vec_size<<" > "<<spacer_len_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG 
        logInfo("\tPassed test 5a. Spacer lengths do not differ too much: "<<abs(sum_spacer_to_spacer_len_diff)/spacer_vec_size<<" < "<<spacer_len_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif    
        int repeat_to_spacer_len_cutoff = spacer_vec_size * CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF;
        
        if (abs(sum_repeat_to_spacer_len_diff) > repeat_to_spacer_len_cutoff) 
        {
//-DDEBUG#ifdef DEBUG
            logInfo("\tFailed test 5b. Repeat to spacer lengths differ too much: "<<abs(sum_repeat_to_spacer_len_diff)/spacer_vec_size<<" > "<<repeat_to_spacer_len_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("\tPassed test 5b. Repeat to spacer lengths do not differ too much: "<<abs(sum_repeat_to_spacer_len_diff)/spacer_vec_size<<" < "<<repeat_to_spacer_len_cutoff/spacer_vec_size, 8);
//-DDEBUG#endif
        
    }
    // short read only one spacer
    else
    {
        std::string spacer = tmp_holder->spacerStringAt(0);
        float similarity = PatternMatcher::getStringSimilarity(repeat, spacer);
        if (similarity > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY) 
        {
//-DDEBUG#ifdef DEBUG
            logInfo("\tFailed test 4. Spacer is too similar to the repeat: "<<similarity<<" > "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("\tPassed test 4. Spacer is not too similar to the repeat: "<<similarity<<" < "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
//-DDEBUG#endif        
        if (abs((int)spacer.length() - (int)repeat.length()) > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF) 
        {
//-DDEBUG#ifdef DEBUG
            logInfo("\tFailed test 5. Repeat to spacer length differ too much: "<<abs((int)spacer.length() - (int)repeat.length())<<" > "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
//-DDEBUG#endif
            return false;
        }
//-DDEBUG#ifdef DEBUG
        logInfo("\tPassed test 5. Repeat to spacer length do not differ too much: "<<abs((int)spacer.length() - (int)repeat.length())<<" < "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
//-DDEBUG#endif
    }
    return true;
    
}

bool isRepeatLowComplexity(std::string& repeat)
{
    int c_count = 0;
    int g_count = 0;
    int a_count = 0;
    int t_count = 0;
    int n_count = 0;
    
    int curr_repeat_length = (int)repeat.length();
    
    int cut_off = (int)(curr_repeat_length * CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD);
    
    std::string::iterator dr_iter = repeat.begin();
    //int i = dr_match.DR_StartPos;
    while (dr_iter != repeat.end()) 
    {
        switch (*dr_iter) 
        {
            case 'c':
            case 'C':
                c_count++; break;
            case 't': 
            case 'T':
                t_count++; break;
            case 'a':
            case 'A':
                a_count++; break;
            case 'g':
            case 'G':
                g_count++; break;
            default: n_count++; break;
        }
        dr_iter++;
    }
    if (a_count > cut_off) return true; 
    else if (t_count > cut_off) return true; 
    else if (g_count > cut_off) return true;
    else if (c_count > cut_off) return true;
    else if (n_count > cut_off) return true;   
    return false;
}

#pragma mark -
#pragma mark ReadHolder Interface

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmpReadholder)
{

    //add the header for the matched readFOG
    
    //logInfo(tmpReadholder->repeatStringAt(2), 8);
    std::string dr_lowlexi = tmpReadholder->DRLowLexi();
    //logInfo(dr_lowlexi, 8);
    StringToken st = mStringCheck->getToken(dr_lowlexi);
    if(0 == st)
    {
        // new guy
        st = mStringCheck->addString(dr_lowlexi);
        (*mReads)[st] = new ReadList();
    }
    (*mReads)[st]->push_back(tmpReadholder);
}

//#pragma mark -
//#pragma mark Utilities
//// turn our map into a vector using just the keys
//void map2Vector(lookupTable& patternsHash, std::vector<std::string>& patterns)
//{
//    
//    lookupTable::iterator iter = patternsHash.begin();
//    while (iter != patternsHash.end()) 
//    {
//        patterns.push_back(iter->first);
//        iter++;
//    }
//}



