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
#include <exception>
#include <ctime>

// local includes
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "WuManber.h"
#include "PatternMatcher.h"
#include "SeqUtils.h"
#include "kseq.h"
#include "StlExt.h"
#include "config.h"
#include "Exception.h"

/* 
declare the type of file handler and the read() function
as described here:
http://lh3lh3.users.sourceforge.net/parsefastq.shtml

THIS JUST DEFINES A BUNCH OF **templated** structs

*/
    KSEQ_INIT(gzFile, gzread)


READ_TYPE decideWhichSearch(const char *inputFastq, float * aveReadLength, const options& opts)
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
        if(read_counter > CRASS_DEF_MAX_READS_FOR_DECISION)
            break;
    }
    kseq_destroy(seq); // destroy seq
    gzclose(fp);
   
    
    *aveReadLength = total_base / read_counter;
    logInfo("Average read length (of the first "<< read_counter<<" reads): "<< *aveReadLength, 2);
    
    // long reads defined by having at least two spacers 4DR + 2SP
    unsigned int long_read_cutoff = (4*opts.lowDRsize) + (2*opts.lowSpacerSize);
    if((total_base / read_counter) > long_read_cutoff)
    {
        return LONG_READ;
    }
    return SHORT_READ;
}


// CRT search
int scanRight(ReadHolder * tmp_holder, std::string& pattern, unsigned int minSpacerLength, unsigned int scanRange)
{
    #ifdef DEBUG
    logInfo("Scanning Right for more repeats:", 9);
    #endif
    unsigned int start_stops_size = tmp_holder->getStartStopListSize();
    
    unsigned int pattern_length = (unsigned int)pattern.length();
    
    // final start index
    unsigned int last_repeat_index = tmp_holder->getRepeatAt(start_stops_size - 2);
    
    //second to final start index
    unsigned int second_last_repeat_index = tmp_holder->getRepeatAt(start_stops_size - 4);
    
    unsigned int repeat_spacing = last_repeat_index - second_last_repeat_index;
    
    #ifdef DEBUG
    logInfo(start_stops_size<<" : "<<pattern_length<<" : "<<last_repeat_index<<" : "<<second_last_repeat_index<<" : "<<repeat_spacing, 9);
    #endif
    
    int candidate_repeat_index, position;
    
    unsigned int begin_search, end_search;
    
    unsigned int read_length = (unsigned int)tmp_holder->getSeqLength();
    bool more_to_search = true;
    while (more_to_search)
    {
        candidate_repeat_index = last_repeat_index + repeat_spacing;
        begin_search = candidate_repeat_index - scanRange;
        end_search = candidate_repeat_index + pattern_length + scanRange;
        #ifdef DEBUG
        logInfo(candidate_repeat_index<<" : "<<begin_search<<" : "<<end_search, 9);
        #endif
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right
        unsigned int scanRightMinBegin = last_repeat_index + pattern_length + minSpacerLength;
        
        if (begin_search < scanRightMinBegin)
        {
            begin_search = scanRightMinBegin;
        }
        if (begin_search > read_length - 1)
        {
            #ifdef DEBUG
            logInfo("returning... "<<begin_search<<" > "<<read_length - 1, 9);
            #endif
            return read_length - 1;
        }
        if (end_search > read_length)
        {
            end_search = read_length;
        }
        
        if ( begin_search >= end_search)
        {
            #ifdef DEBUG
            logInfo("Returning... "<<begin_search<<" >= "<<end_search, 9);
            #endif
            return end_search;
        }
        /******************** end range checks ********************/
        
        std::string text = tmp_holder->substr(begin_search, (end_search - begin_search));
        
        #ifdef DEBUG
        logInfo(pattern<<" : "<<text, 9);
        #endif
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

int longReadSearch(const char * inputFastq, const options& opts, ReadMap * mReads, StringCheck * mStringCheck, lookupTable& patternsHash, lookupTable& readsFound)
{
    //-----
    // Code lifted from CRT, ported by Connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, log_counter = 0;
    static int read_counter = 0;
    // initialize seq
    seq = kseq_init(fp);
    //ReadHolder * tmp_holder = new ReadHolder();
    // read sequence 
    time_t time_start, time_current;
    time(&time_start);
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        if (log_counter == CRASS_DEF_READ_COUNTER_LOGGER) 
        {
            time(&time_current);
            double diff = difftime(time_current, time_start);
            time_start = time_current;
            std::cout<<"["<<PACKAGE_NAME<<"_longReadFinder]: "<< "Processed "<<read_counter<<" ...";
            std::cout<<diff<<"sec"<<std::endl;
            log_counter = 0;
        }
        
        // grab a readholder
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
        // test if it has a comment entry and a quality entry
        if (seq->comment.s) 
        {
            tmp_holder->setComment(seq->comment.s);
        }
        if (seq->qual.s) 
        {
            tmp_holder->setQual(seq->qual.s);
        }
        
        bool match_found = false;

        if (opts.removeHomopolymers)
        {
            // RLE is necessary...
            tmp_holder->encode();
        }
        std::string read = tmp_holder->getSeq();
        
        // get the length of this sequence
        unsigned int seq_length = (unsigned int)read.length();
        

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
            delete tmp_holder;
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
            
            std::string text;
            std::string pattern;
            try {
                text = read.substr(beginSearch, (endSearch - beginSearch));

            } catch (std::exception& e) {
                throw crispr::substring_exception(e.what(), text.c_str(), beginSearch, (endSearch - beginSearch), __FILE__, __LINE__, __PRETTY_FUNCTION__);
            }
            try {
                pattern = read.substr(j, opts.searchWindowLength);
            } catch (std::exception& e) {
                throw crispr::substring_exception(e.what(), read.c_str(), j, opts.searchWindowLength, __FILE__, __LINE__, __PRETTY_FUNCTION__);
            }

            //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
            int pattern_in_text_index = -1;
                pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);

            if (pattern_in_text_index >= 0)
            {
                tmp_holder->startStopsAdd(j,  j + opts.searchWindowLength);
                unsigned int found_pattern_start_index = beginSearch + (unsigned int)pattern_in_text_index;
                
                tmp_holder->startStopsAdd(found_pattern_start_index, found_pattern_start_index + opts.searchWindowLength);
                scanRight(tmp_holder, pattern, opts.lowSpacerSize, 24);
            }

            if ( (tmp_holder->numRepeats() > opts.minNumRepeats) ) //tmp_holder->numRepeats is half the size of the StartStopList
            {
#ifdef DEBUG
                logInfo(tmp_holder->getHeader(), 8);
                logInfo("\tPassed test 1. More than "<<opts.minNumRepeats<< " ("<<tmp_holder->numRepeats()<<") repeated kmers found", 8);
#endif

                unsigned int actual_repeat_length = extendPreRepeat(tmp_holder, opts.searchWindowLength);

                if ( (actual_repeat_length >= opts.lowDRsize) && (actual_repeat_length <= opts.highDRsize) )
                {
#ifdef DEBUG

                    logInfo("\tPassed test 2. The repeat length is between "<<opts.lowDRsize<<" and "<<opts.highDRsize, 8);
#endif
                    if (opts.removeHomopolymers) 
                    {
                        tmp_holder->decode();
                    }
                    
                    // drop partials
                    tmp_holder->dropPartials();
                    if (qcFoundRepeats(tmp_holder, opts.lowSpacerSize, opts.highSpacerSize))
                    {
#ifdef DEBUG
                        logInfo("Passed all tests!", 8);
                        logInfo("Potential CRISPR containing read found: "<<tmp_holder->getHeader(), 7);
                        logInfo(tmp_holder->getSeq(), 9);
                        logInfo("-------------------", 7)
#endif                            

                        //ReadHolder * candidate_read = new ReadHolder();
                        //*candidate_read = *tmp_holder;
                        //addReadHolder(mReads, mStringCheck, candidate_read);
                        addReadHolder(mReads, mStringCheck, tmp_holder);
                        match_found = true;
                        patternsHash[tmp_holder->repeatStringAt(0)] = true;
                        readsFound[tmp_holder->getHeader()] = true;
                        break;
                    }
                }
#ifdef DEBUG                
                else
                {
                    logInfo("\tFailed test 2. Repeat length: "<<tmp_holder->getRepeatLength() << " : " << match_found, 8); 
                }
#endif
                j = tmp_holder->back() - 1;
            }
            tmp_holder->clearStartStops();
        }
        if (!match_found) 
        {
            delete tmp_holder;
        }
        log_counter++;
        read_counter++;
    }
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    //delete tmp_holder;
    logInfo("finished processing file:"<<inputFastq, 1);    
    time(&time_current);
    double diff = difftime(time_current, time_start);
    time_start = time_current;
    std::cout<<"["<<PACKAGE_NAME<<"_longReadFinder]: "<< "Processed "<<read_counter<<" ...";
    std::cout<<diff<<std::endl;
    logInfo("So far " << mReads->size()<<" direct repeat variants have been found from " << read_counter << " reads", 2);
    return 0;
}

int shortReadSearch(const char * inputFastq, const options& opts, lookupTable& patternsHash, lookupTable& readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, log_counter = 0;
    static int read_counter = 0;
    // initialize seq
    seq = kseq_init(fp);
    time_t time_start, time_current;
    time(&time_start);
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        if (log_counter == CRASS_DEF_READ_COUNTER_LOGGER) 
        {
            time(&time_current);
            double diff = difftime(time_current,time_start );
            time_start = time_current;
            std::cout<<"["<<PACKAGE_NAME<<"_shortReadFinder]: "<<"Processed "<<read_counter<<" ...";
            //std::cout.setf(std::ios::fixed);
            //std::cout.precision(2);
            std::cout<<diff<<"sec"<<std::endl;
            log_counter = 0;
        }
        
        bool match_found = false;

        // create the read holder
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
        if (seq->comment.s) 
        {
            tmp_holder->setComment(seq->comment.s);
        }
        if (seq->qual.s) 
        {
            tmp_holder->setQual(seq->qual.s);
        }
        
        if (opts.removeHomopolymers)
        {
            // RLE is necessary...
            tmp_holder->encode();
        }
        
        std::string read = tmp_holder->getSeq();

        unsigned int seq_length = (unsigned int)read.length();
        unsigned int search_end = seq_length - opts.lowDRsize - 1;
        unsigned int final_index = seq_length - 1;
        
        
        // boyer-moore search
        for (unsigned int first_start = 0; first_start < search_end; first_start++)
        {
            unsigned int search_begin = first_start + opts.lowDRsize + opts.lowSpacerSize;
            
            if (search_begin >= search_end ) break;
            
            // do the search
            int second_start = -1;

            try {
                second_start = PatternMatcher::bmpSearch( read.substr(search_begin), read.substr(first_start, opts.lowDRsize) );
            } catch (std::exception& e) {
                throw crispr::exception( __FILE__, __LINE__, __PRETTY_FUNCTION__,e.what());
            }

            // check to see if we found something
            if (second_start > -1) 
            {
                // bingo!
                second_start += search_begin;
                unsigned int second_end = (unsigned int)second_start + opts.lowDRsize;
                unsigned int first_end = first_start + opts.lowDRsize;
  
                unsigned int next_index = second_end + 1;
                
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                if ( next_index <= final_index) 
                {
                    // read through the subsuquent bases untill they don't match
                    unsigned int extenstion_length = 0;
                    while (read.at(first_end + extenstion_length) == read.at(second_end + extenstion_length)) 
                    {
#ifdef DEBUG
                        logInfo(read.at(first_end + extenstion_length)<<" == "<<read.at(second_end + extenstion_length),20);
#endif
                        extenstion_length++;
                        next_index++;
                        if (next_index > final_index) break;

                    }
#ifdef DEBUG
                    logInfo("A: FS: " << first_start<< " FE: "<<first_end<<" SS: "<< second_start<< " SE: "<<second_end<<" EX: "<<extenstion_length, 8);
#endif
                    tmp_holder->startStopsAdd(first_start, first_end + extenstion_length);
                    tmp_holder->startStopsAdd(second_start, second_end + extenstion_length);
                    tmp_holder->setRepeatLength(second_end - second_start + extenstion_length);
                }
                else
                {
#ifdef DEBUG
                    logInfo("B: FS: " << first_start<< " FE: "<<first_end<<" SS: "<< second_start<< " SE: "<<second_end, 8);
#endif
                    tmp_holder->startStopsAdd(first_start, first_end);
                    tmp_holder->startStopsAdd(second_start, second_end);
                    tmp_holder->setRepeatLength(second_end - second_start);
                }
                
                // the low side will always be true since we search for the lowDRSize
                if ( tmp_holder->getRepeatLength() <= opts.highDRsize )
                {

                    if ((tmp_holder->getAverageSpacerLength() >= opts.lowSpacerSize) && (tmp_holder->getAverageSpacerLength() <= opts.highSpacerSize))
                    {
                        if (qcFoundRepeats(tmp_holder, opts.lowSpacerSize, opts.highSpacerSize)) 
                        {

                            match_found = true;
#ifdef DEBUG
                            logInfo("Potential CRISPR containing read found: "<<tmp_holder->getHeader(), 7);
                            logInfo(read, 9);
                            logInfo("-------------------", 7)
#endif
                            patternsHash[tmp_holder->repeatStringAt(0)] = true;
                            readsFound[tmp_holder->getHeader()] = true;
                            addReadHolder(mReads, mStringCheck, tmp_holder);
                            break;
                        }
                    }
#ifdef DEBUG
                    else
                    {
                        logInfo("\tFailed test 2. The spacer length is not between "<<opts.lowSpacerSize<<" and "<<opts.highSpacerSize<<": "<<tmp_holder->getAverageSpacerLength(), 8);
                    }
#endif
                }
#ifdef DEBUG
                else
                {
                    logInfo("\tFailed test 1. The repeat length is larger than "<<opts.highDRsize<<": " << tmp_holder->getRepeatLength(), 8);
                }
#endif
                first_start = tmp_holder->back();
            }
        }
        if (!match_found)
        {
            delete tmp_holder;
        }
        log_counter++;
        read_counter++;
    }

    // clean up
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("Finished processing file:"<<inputFastq, 1);
    time(&time_current);
    double diff = difftime( time_current, time_start);
    time_start = time_current;
    std::cout<<"["<<PACKAGE_NAME<<"_shortReadFinder]: "<<"Processed "<<read_counter<<" ...";
    std::cout<<diff<<std::endl;    logInfo("So far " << mReads->size()<<" direct repeat variants have been found", 2);
    return 0;
}


void findSingletons(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    //-----
    // given a potentially large set of patterns, call an innefficent function
    // on a possible broken down subset
    //
    std::vector<std::string> patterns;
    mapToVector(patternsHash, patterns);
    try
    {
        if (patterns.empty())
        {
            logError("No patterns in vector for multimatch");
            throw "No patterns in vector for multimatch";
        }
    }
    catch (char * c)
    {
        std::cerr << c << std::endl;
        return;
    }   

    // If the patterns vector is too long, we'll need to break it up!
    // in any case, we need to pass through a vector of vectors
    std::vector<std::vector<std::string> * > vec_vec_patterns;
    int cut_start = 0;
    int cut_length = CRASS_DEF_MAX_SING_PATTERNS;
    int total_pattern_size = (int)patterns.size();
    std::vector<std::string>::iterator start_iter, end_iter;
    while(total_pattern_size > 0)
    {
        if(total_pattern_size <  cut_length)
        {
            cut_length = total_pattern_size;
        }
        
        // in the second last iteration, we should check to
        // see if halving is a better option
        if((total_pattern_size < 2 * CRASS_DEF_MAX_SING_PATTERNS) && (total_pattern_size > CRASS_DEF_MAX_SING_PATTERNS))
        {
            // in most cases, halving is a better option
            cut_length = (int)(total_pattern_size/2);
        }
        
        // make out iterators
        start_iter = patterns.begin() + cut_start;
        end_iter = start_iter + cut_length;
        
        // get a new vector to work with
        // intialised with the subset of the patterns
        std::vector<std::string> * tmp_vec = new std::vector<std::string>();
        tmp_vec->insert(tmp_vec->begin(), start_iter, end_iter);
        vec_vec_patterns.push_back(tmp_vec);
        
        // update our counters
        cut_start += cut_length;
        total_pattern_size -= cut_length;
        cut_length = CRASS_DEF_MAX_SING_PATTERNS;
    }

    findSingletonsMultiVector(inputFastq, opts, vec_vec_patterns, readsFound, mReads, mStringCheck);
    
    // clean up
    std::vector<std::vector<std::string> * >::iterator vv_iter = vec_vec_patterns.begin();
    std::vector<std::vector<std::string> * >::iterator vv_last = vec_vec_patterns.end();
    while(vv_iter != vv_last)
    {
        if(*vv_iter != NULL)
        {
            delete *vv_iter;
            *vv_iter = NULL;
        }
        vv_iter++;
    }
    
}

void findSingletonsMultiVector(const char *inputFastq, const options &opts, std::vector<std::vector<std::string> *> &patterns, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    //-----
    // Find sings given a vector of vectors of patterns
    //


    // make a wumanber for each set of patterns
    std::vector<WuManber*> wu_mans;
    std::vector<std::vector<std::string> *>::iterator pats_iter = patterns.begin();
    std::vector<std::vector<std::string> *>::iterator pats_last = patterns.end();
    while(pats_iter != pats_last)
    {
        WuManber * search = new WuManber();
        search->Initialize(*(*pats_iter));
        wu_mans.push_back(search);
        pats_iter++;
    }
    
    
    // now we got lots of wumanbers, search each string

    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    
    seq = kseq_init(fp);

    int l;
    int log_counter = 0;
    static int read_counter = 0;
    int old_number = (int)mReads->size();

    std::vector<WuManber*>::iterator wm_iter =  wu_mans.begin();
    std::vector<WuManber*>::iterator wm_last =  wu_mans.end();
    time_t time_start, time_current;
    time(&time_start);
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        // seq is a read what we love
        // search it for the patterns until found
        bool found = false;
        if (log_counter == CRASS_DEF_READ_COUNTER_LOGGER) 
        {
            time(&time_current);
            double diff = difftime(time_current, time_start);
            time_start = time_current;
            std::cout<<"["<<PACKAGE_NAME<<"_singletonFinder]: "<<"Processed "<<read_counter<<" ...";
            //std::cout.setf(std::ios::fixed);
            //std::cout.precision(2);
            std::cout<<diff<<"sec"<<std::endl;
            log_counter = 0;
        }
        // reset these mofos
        wm_iter =  wu_mans.begin();
        wm_last =  wu_mans.end();
        pats_iter = patterns.begin();

        while(wm_iter != wm_last)
        {
            ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
            if (seq->comment.s) 
            {
                tmp_holder->setComment(seq->comment.s);
            }
            if (seq->qual.s) 
            {
                tmp_holder->setQual(seq->qual.s);
            }
            
            if (opts.removeHomopolymers) 
            {
                tmp_holder->encode();
            }
            std::string read = tmp_holder->getSeq();
            
            //initialize with an impossible number
            int start_pos = -1;
            std::string found_repeat = (*wm_iter)->Search(read.length(), read.c_str(), *(*pats_iter), start_pos);
            
            if (start_pos != -1)
            {
                if (readsFound.find(tmp_holder->getHeader()) == readsFound.end())
                {
    #ifdef DEBUG
                    logInfo("new read recruited: "<<tmp_holder->getHeader(), 9);
                    logInfo(tmp_holder->getSeq(), 10);
    #endif
                    unsigned int DR_end = (unsigned int)start_pos + (unsigned int)found_repeat.length();
                    if(DR_end >= (unsigned int)read.length())
                    {
                        DR_end = (unsigned int)read.length() - 1;
                    }
                    tmp_holder->startStopsAdd(start_pos, DR_end);
                    addReadHolder(mReads, mStringCheck, tmp_holder);
                    found = true;
                }
                else
                    delete tmp_holder;
            }
            else
            {
                delete tmp_holder;
            }
            
            if(found)
                break;
            
            pats_iter++;
            wm_iter++;
        }
        log_counter++;
        read_counter++;
    }
    
    // clean up
    wm_iter =  wu_mans.begin();
    wm_last =  wu_mans.end();
    while(wm_iter != wm_last)
    {
        if(*wm_iter != NULL)
        {
            delete *wm_iter;
            *wm_iter = NULL;
        }
        wm_iter++;
    }
    time(&time_current);
    double diff = difftime(time_current, time_start);
    time_start = time_current;
    std::cout<<"["<<PACKAGE_NAME<<"_singletonFinder]: "<<"Processed "<<read_counter<<" ...";
    std::cout<<diff<<std::endl;
    logInfo("Finished second iteration. An extra " << mReads->size() - old_number<<" variants were recruited", 2);
}

unsigned int extendPreRepeat(ReadHolder * tmp_holder, int searchWindowLength)
{
#ifdef DEBUG
    logInfo("Extending Prerepeat...", 9);
#endif
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
#ifdef DEBUG
    logInfo("cutoff: "<<cut_off, 9);
#endif
    
    
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
#ifdef DEBUG
        logInfo(i<<" : "<<curr_repeat_spacing, 10);
#endif
        
        // if it is shorter than what we already have, make it the shortest
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }
#ifdef DEBUG
    logInfo("shortest repeat spacing: "<<shortest_repeat_spacing, 9);
#endif
    unsigned int right_extension_length = 0;
    // don't search too far  
    unsigned int max_right_extension_length = shortest_repeat_spacing;// - minSpacerLength;
#ifdef DEBUG
    logInfo("max right extension length: "<<max_right_extension_length, 9);
#endif
    // Sometimes we shouldn't use the far right DR. (it may lie too close to the end)
    unsigned int DR_index_end = end_index;
    unsigned int dist_to_end = tmp_holder->getSeqLength() - last_repeat_start_index - 1;
    if(dist_to_end < max_right_extension_length)
    {
#ifdef DEBUG
        logInfo("removing end partial: "<<dist_to_end<<" < "<<max_right_extension_length, 9);
#endif
        DR_index_end -= 2;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off) 
        {
            cut_off =2;
        }
#ifdef DEBUG
        logInfo("new cutoff: "<<cut_off, 9);
#endif
    }
    std::string curr_repeat;
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
    while (max_right_extension_length > 0)
    {
        for (unsigned int k = 0; k < DR_index_end; k+=2 )
        {
#ifdef DEBUG
            logInfo(k<<" : "<<tmp_holder->getRepeatAt(k) + tmp_holder->getRepeatLength(), 10);
#endif
            // look at the character just past the end of the last repeat
            // make sure our indicies make some sense!
            if((tmp_holder->getRepeatAt(k) + tmp_holder->getRepeatLength()) >= (unsigned int)tmp_holder->getSeqLength())
            {    
                k = DR_index_end;
            }
            else
            {
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
        }
        
#ifdef DEBUG
        logInfo("R: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << tmp_holder->getRepeatLength() << " : " << max_right_extension_length, 9);
#endif
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
    int test_for_negative = shortest_repeat_spacing - tmp_holder->getRepeatLength();// - minSpacerLength;
    unsigned int max_left_extension_length = (test_for_negative >= 0)? (unsigned int) test_for_negative : 0;

#ifdef DEBUG
    logInfo("max left extension length: "<<max_left_extension_length, 9);
#endif
    // and again, we may not wat to use the first DR
    unsigned int DR_index_start = 0;
    if(max_left_extension_length > first_repeat_start_index)
    {
#ifdef DEBUG
        logInfo("removing start partial: "<<max_left_extension_length<<" > "<<first_repeat_start_index, 9);
#endif
        DR_index_start+=2;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off) 
        {
            cut_off = 2;
        }
#ifdef DEBUG
        logInfo("new cutoff: "<<cut_off, 9);
#endif
    }
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (left_extension_length < max_left_extension_length)
    {
        for (unsigned int k = DR_index_start; k < end_index; k+=2 )
        {
#ifdef DEBUG
            logInfo(k<<" : "<<tmp_holder->getRepeatAt(k) - left_extension_length - 1, 10);
#endif
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
#ifdef DEBUG
        logInfo("L:" << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << tmp_holder->getRepeatLength() << " : " << left_extension_length, 9);
#endif
        
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
#ifdef DEBUG    
    logInfo("Repeat positions:", 9);
#endif
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
#ifdef DEBUG    
        logInfo(*repeat_iter<<","<<*(repeat_iter+1), 9);
#endif
        repeat_iter += 2;
    }
    
    return (unsigned int)tmp_holder->getRepeatLength();
    
}


//need at least two elements
bool qcFoundRepeats(ReadHolder * tmp_holder, int minSpacerLength, int maxSpacerLength)
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
#ifdef DEBUG
        logInfo("\tFailed test 3. The repeat is low complexity", 8);
#endif
        return false;
    }
    
#ifdef DEBUG
    logInfo("\tPassed test 3. The repeat is not low complexity", 8);
#endif
    // test for a long or short read
    int single_compare_index = 0;
    bool is_short = (2 > tmp_holder->numSpacers()); 
    if(!is_short) 
    {
        // for holding stats
        float ave_spacer_to_spacer_len_difference = 0.0;
        float ave_repeat_to_spacer_len_difference = 0.0;
        float ave_spacer_to_spacer_difference = 0.0;
        float ave_repeat_to_spacer_difference = 0.0;
        int min_spacer_length = 10000000;
        int max_spacer_length = 0;
        int num_compared = 0;
        
        // we need to do things a little differently depending on whether of not this guy starts 
        // or ends on a spacer...
        int sp_start_offset = 0;
        int sp_end_offset = 1;
        if (tmp_holder->startStopsAt(0) != 0)
        {
            // starts with a spacer
            //MI std::cout << "sp_start" << std::endl;
            sp_start_offset++;
        }
        if (tmp_holder->getSeqLength() != (int)tmp_holder->back() + 1) 
        {
            // ends with a spacer
            //MI std::cout << "sp_end" << std::endl;
            sp_end_offset++;
        }

        // now go through the spacers and check for similarities
        std::vector<std::string> spacer_vec = tmp_holder->getAllSpacerStrings();
        std::vector<std::string>::iterator spacer_iter = spacer_vec.begin();
        std::vector<std::string>::iterator spacer_last = spacer_vec.end();
        spacer_iter += sp_start_offset;
        spacer_last -= sp_end_offset;
        while (spacer_iter != spacer_last) 
        {
            num_compared++;
            ave_repeat_to_spacer_difference += PatternMatcher::getStringSimilarity(repeat, *spacer_iter);
            float ss_diff = PatternMatcher::getStringSimilarity(*spacer_iter, *(spacer_iter + 1));
            ave_spacer_to_spacer_difference += ss_diff;
            //MI std::cout << ss_diff << " : " << *spacer_iter << " : " << *(spacer_iter + 1) << std::endl;
            ave_spacer_to_spacer_len_difference += ((float)(spacer_iter->size()) - (float)((spacer_iter + 1)->size()));
            ave_repeat_to_spacer_len_difference +=  ((float)(repeat.size()) - (float)(spacer_iter->size()));
            spacer_iter++;
        }

        // now look for max and min lengths!
        spacer_iter = spacer_vec.begin() + sp_start_offset;
        spacer_last++;
        while (spacer_iter != spacer_last) 
        {
            num_compared++;
            if((int)(spacer_iter->length()) < min_spacer_length)
            {
                min_spacer_length = (int)(spacer_iter->length()); 
            }
            if((int)(spacer_iter->length()) > max_spacer_length)
            {
                max_spacer_length = (int)(spacer_iter->length()); 
            }
            spacer_iter++;
        }

        // we may not have compared anything...
        if(num_compared == 0)
        {
            // there are only three spacers in this read and the read begins and ends on a spacer
            // we will still need to do some comparisons
            is_short = true;
            single_compare_index = 1;
        }
        else
        {
            // we can keep going!
            ave_spacer_to_spacer_difference /= num_compared;
            ave_repeat_to_spacer_difference /= num_compared;
            ave_spacer_to_spacer_len_difference = abs(ave_spacer_to_spacer_len_difference /= num_compared);
            ave_repeat_to_spacer_len_difference = abs(ave_repeat_to_spacer_len_difference /= num_compared);
            
            /*
             * MAX AND MIN SPACER LENGTHS
             */
            if (min_spacer_length < minSpacerLength) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 4a. Min spacer length out of range: "<<min_spacer_length<<" < "<<minSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4a. Min spacer length within range: "<<min_spacer_length<<" > "<<minSpacerLength, 8);
#endif
            if (min_spacer_length > maxSpacerLength) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 4b. Max spacer length out of range: "<<max_spacer_length<<" > "<<maxSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4b. Max spacer length within range: "<<max_spacer_length<<" < "<<maxSpacerLength, 8);
#endif
            
            /*
             * REPEAT AND SPACER CONTENT SIMILARITIES
             */
            if (ave_spacer_to_spacer_difference > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 5a. Spacers are too similar: "<<ave_spacer_to_spacer_difference<<" > "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 5a. Spacers are not too similar: "<<ave_spacer_to_spacer_difference<<" < "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif    
            if (ave_repeat_to_spacer_difference > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY)
            {
#ifdef DEBUG
                logInfo("\tFailed test 5b. Spacers are too similar to the repeat: "<<ave_repeat_to_spacer_difference<<" > "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 5b. Spacers are not too similar to the repeat: "<<ave_repeat_to_spacer_difference<<" < "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif    
            
            /*
             * REPEAT AND SPACER LENGTH SIMILARITIES
             */
            if (ave_spacer_to_spacer_len_difference > CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF) 
            {
#ifdef DEBUG 
                logInfo("\tFailed test 6a. Spacer lengths differ too much: "<<ave_spacer_to_spacer_len_difference<<" > "<<CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF, 8);
#endif
                return false;
            }
#ifdef DEBUG 
            logInfo("\tPassed test 6a. Spacer lengths do not differ too much: "<<ave_spacer_to_spacer_len_difference<<" < "<<CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF, 8);
#endif    
            if (ave_repeat_to_spacer_len_difference > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 6b. Repeat to spacer lengths differ too much: "<<ave_repeat_to_spacer_len_difference<<" > "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 6b. Repeat to spacer lengths do not differ too much: "<<ave_repeat_to_spacer_len_difference<<" < "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
#endif
            
        }
    }
    
    // Are we testing a short read or only one spacer?
    if(is_short)
    {
        std::string spacer = tmp_holder->spacerStringAt(single_compare_index);
        float similarity = PatternMatcher::getStringSimilarity(repeat, spacer);
        if (similarity > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY) 
        {
            /*
             * MAX AND MIN SPACER LENGTHS
             */
            if ((int)spacer.length() < minSpacerLength) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 4a. Min spacer length out of range: "<<spacer.length()<<" < "<<minSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4a. Min spacer length within range: "<<spacer.length()<<" > "<<minSpacerLength, 8);
#endif
            if ((int)spacer.length() > maxSpacerLength) 
            {
#ifdef DEBUG
                logInfo("\tFailed test 4b. Max spacer length out of range: "<<spacer.length()<<" > "<<maxSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
             logInfo("\tPassed test 4b. Max spacer length within range: "<<spacer.length()<<" < "<<maxSpacerLength, 8);
#endif
            /*
             * REPEAT AND SPACER CONTENT SIMILARITIES
             */ 
#ifdef DEBUG
            logInfo("\tFailed test 5. Spacer is too similar to the repeat: "<<similarity<<" > "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif
            return false;
        }
#ifdef DEBUG
        logInfo("\tPassed test 5. Spacer is not too similar to the repeat: "<<similarity<<" < "<<CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
#endif        
        
        /*
         * REPEAT AND SPACER LENGTH SIMILARITIES
         */
        if (abs((int)spacer.length() - (int)repeat.length()) > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF) 
        {
#ifdef DEBUG
            logInfo("\tFailed test 6. Repeat to spacer length differ too much: "<<abs((int)spacer.length() - (int)repeat.length())<<" > "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
#endif
            return false;
        }
#ifdef DEBUG
        logInfo("\tPassed test 6. Repeat to spacer length do not differ too much: "<<abs((int)spacer.length() - (int)repeat.length())<<" < "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
#endif
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




void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmpReadholder)
{


    std::string dr_lowlexi = tmpReadholder->DRLowLexi();
#ifdef DEBUG
    logInfo("Direct repeat: "<<dr_lowlexi, 10);
#endif
    StringToken st = mStringCheck->getToken(dr_lowlexi);

    if(0 == st)
    {
        // new guy
        st = mStringCheck->addString(dr_lowlexi);
        (*mReads)[st] = new ReadList();
    }
#ifdef DEBUG
    logInfo("String Token: "<<st, 10);
#endif
    (*mReads)[st]->push_back(tmpReadholder);
}

