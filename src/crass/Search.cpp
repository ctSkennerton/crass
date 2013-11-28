// File: Search.cpp
// Original Author: Connor Skennerton 2013
// Based on libcrispr.cpp by Mike Imelfort
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Searches for CRISPRs in reads or Files
//
// --------------------------------------------------------------------
//  Copyright  2011-2013 Michael Imelfort and Connor Skennerton
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
#include <map>
// local includes
#include "Search.h"
#include "PatternMatcher.h"

int crass::Search::searchFileSerial(const char *inputFastq)

{
    //-----
    // Wrapper used for searching reads for DRs
    // depending on the length of the read.
	// this funciton may use the boyer moore algorithm
    // or the CRT search algorithm
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t * seq = kseq_init(fp);
    
    int l, log_counter, max_read_length;
    log_counter = max_read_length = 0;
    static int read_counter = 0;
    time_t time_current;
    
    int length_cutoff = (2 * mMinRepeatLength) + mMinSpacerLength;
    int reads_found_with_crispr = 0;
    while ( (l = kseq_read(seq)) >= 0 )
    {
        max_read_length = (l > max_read_length) ? l : max_read_length;
        /*if (log_counter == CRASS_DEF_READ_COUNTER_LOGGER)
        {
            time(&time_current);
            double diff = difftime(time_current, time_start);
            //time_start = time_current;
            
            std::cout<<"\r["<<PACKAGE_NAME<<"_patternFinder]: "
            << "Processed "<<read_counter<<" ...";
            std::cout<<diff<<" sec"<<std::flush;
            log_counter = 0;
        }*/
        crass::RawRead read(seq->name.s, "", seq->seq.s, "", crass::RepeatArray());
        if (seq->comment.s != NULL) {
            read.comment(seq->comment.s);
        }
        if(seq->qual.s != NULL) {
            read.quality(seq->qual.s);
        }
#if SEARCH_SINGLETON
            SearchCheckerList::iterator debug_iter = debugger->find(seq->name.s);
            if (debug_iter != debugger->end()) {
                changeLogLevel(10);
                std::cout<<"Processing interesting read: "<<debug_iter->first<<std::endl;
            } else {
                changeLogLevel(opts.logLevel);
            }
#endif
            
            if (l > length_cutoff) {
                // perform long read seqrch
                if(readSearch(read))
                    reads_found_with_crispr++;
            }

        log_counter++;
        read_counter++;
    }
    
    kseq_destroy(seq); // destroy seq
    gzclose(fp);
    
    /*logInfo("finished processing file:"<<inputFastq, 1);
    time(&time_current);
    double diff = difftime(time_current, time_start);
    //time_start = time_current;
    std::cout<<"\r["<<PACKAGE_NAME<<"_patternFinder]: "<< "Processed "<<read_counter<<" ...";
    std::cout<<diff<<" sec"<<std::flush;
    logInfo("So far " << mReads->size()<<" direct repeat variants have been found from " << read_counter << " reads", 2);
    */
    logInfo("reads found with crispr: "<<reads_found_with_crispr<<" ("<<(float)reads_found_with_crispr/(float)read_counter<<")" , 1);
    return max_read_length;
}

bool crass::Search::readSearch(crass::RawRead& read)
{

    int seq_length = static_cast<int>(read.seq().length());
    
    
    //the number of bases that can be skipped while we still guarantee that the entire search
    //window will at some point in its iteration thru the sequence will not miss a any repeat
    unsigned int skips = mMinRepeatLength - (2 * mSearchWindowLength - 1);
    if (skips < 1)
    {
        skips = 1;
    }
    
    int searchEnd = seq_length - mMinRepeatLength - mMinSpacerLength - mSearchWindowLength;
    
    if (searchEnd < 0)
    {
        logInfo("Read:"<<read.identifier()<<" length is less than " << mMinRepeatLength + mMinSpacerLength + mSearchWindowLength<<"bp", 1);
        return false;
    }
    
    for ( int j = 0; j <= searchEnd; j = j + skips)
    {
        
        int beginSearch = j + mMinRepeatLength + mMinSpacerLength;
        int endSearch = j + mMaxRepeatLength + mMaxSpacerLength + mSearchWindowLength;


        if (endSearch >= seq_length)
        {
            endSearch = seq_length;
        }
        //should never occur
        if (endSearch < beginSearch)
        {
            endSearch = beginSearch;
        }
        std::string text;
        std::string pattern;
        text = read.seq().substr(beginSearch, (endSearch - beginSearch));

        pattern = read.seq().substr(j, mSearchWindowLength);

        
        int pattern_in_text_index = -1;
        pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);
        
        if (pattern_in_text_index >= 0)
        {
            read.push_back(j,  j + mSearchWindowLength);
            int found_pattern_start_index = beginSearch + pattern_in_text_index;
            
            read.push_back(found_pattern_start_index, found_pattern_start_index + mSearchWindowLength);
            scanRight(read, pattern);
        }
        
        if ( (read.numberOfRepeats() >= mMinSeedCount) )
        {
#ifdef DEBUG
            logInfo(read.identifier(), 8);
            logInfo("\tPassed test 1. At least "<<mMinSeedCount<< " ("<<read.numberOfRepeats()<<") repeated kmers found", 8);
#endif
            
            int actual_repeat_length = extendPreRepeat(read);
            
            if ( (actual_repeat_length >= mMinRepeatLength) && (actual_repeat_length <= mMaxRepeatLength) )
            {
#ifdef DEBUG
                
                logInfo("\tPassed test 2. The repeat length is "<<mMinRepeatLength<<" >= "<< actual_repeat_length <<" <= "<<mMaxRepeatLength, 8);
#endif
                
                if (qcFoundRepeats(read))
                {
#ifdef DEBUG
                    logInfo("Passed all tests!", 8);
                    read.inspect();
                    logInfo("-------------------", 8)
#endif
                    
                    return true;
                }
            }
#ifdef DEBUG
            else
            {
                logInfo("\tFailed test 2. Repeat length: "<<actual_repeat_length/* << " : " << match_found*/, 8);
            }
#endif
            j = (*(read.repeatBegin())).second;
        }
        read.clearRepeatPositions();
    }
    return false;
}

bool crass::Search::readPairSearch(crass::RawRead& read1, crass::RawRead& read2)
{

    if (!readPairSearchCore(read1, read2)) {
        read1.clearRepeatPositions();
        read2.clearRepeatPositions();
        
        if (!readPairSearchCore(read2, read1)) {
            read1.clearRepeatPositions();
            read2.clearRepeatPositions();
            return false;
        } 
    } 
    return true;
}


bool crass::Search::readPairSearchCore(crass::RawRead &read1, crass::RawRead &read2) {
    
    int seq_length = static_cast<int>(read1.seq().length());
    
    
    //the number of bases that can be skipped while we still guarantee that the entire search
    //window will at some point in its iteration thru the sequence will not miss a any repeat
    unsigned int skips = mMinRepeatLength - (2 * mSearchWindowLength - 1);
    if (skips < 1)
    {
        skips = 1;
    }
    
    int searchEnd = seq_length - mMinRepeatLength - mMinSpacerLength - mSearchWindowLength;
    
    if (searchEnd < 0)
    {
        logInfo("Read:"<<read1.identifier()<<" length is less than " << mMinRepeatLength + mMinSpacerLength + mSearchWindowLength<<"bp", 1);
        return false;
    }
    
    for ( int j = 0; j <= searchEnd; j = j + skips)
    {
        
        int beginSearch = j + mMinRepeatLength + mMinSpacerLength;
        int endSearch = j + mMaxRepeatLength + mMaxSpacerLength + mSearchWindowLength;
        
        if (endSearch >= seq_length)
        {
            endSearch = seq_length;
        }
        //should never occur
        if (endSearch < beginSearch)
        {
            endSearch = beginSearch;
        }
        
        std::string text;
        std::string pattern;
        text = read1.seq().substr(beginSearch, (endSearch - beginSearch));
        
        pattern = read1.seq().substr(j, mSearchWindowLength);
        
        
        int pattern_in_text_index = -1;
        pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);
        
        if (pattern_in_text_index >= 0)
        {
            read1.push_back(j,  j + mSearchWindowLength);
            int found_pattern_start_index = beginSearch + pattern_in_text_index;
            
            read1.push_back(found_pattern_start_index, found_pattern_start_index + mSearchWindowLength);
            scanRight(read1, pattern);
            
        }
        
        if ( (read1.numberOfRepeats() >= mMinSeedCount) )
        {
#ifdef DEBUG
            logInfo(read1.identifier(), 8);
            logInfo("\tPassed test 1. At least "<<mMinSeedCount<< " ("<<read1.numberOfRepeats()<<") repeated kmers found", 8);
#endif
            // find seed in second read
            std::vector<int> sposs;
            PatternMatcher::bmpMultiSearch(read2.seq(), pattern, sposs);
            int actual_repeat_length;
            if (!sposs.empty()) {
                for(auto it = sposs.begin(); it != sposs.end(); ++it) {
                    read2.push_back(*it, (*it) + mSearchWindowLength);
                }
                actual_repeat_length = extendPreRepeatPair(read1, read2);
            } else {
                actual_repeat_length = extendPreRepeat(read1);
            }
            
            
            if ( (actual_repeat_length >= mMinRepeatLength) && (actual_repeat_length <= mMaxRepeatLength) )
            {
#ifdef DEBUG
                
                logInfo("\tPassed test 2. The repeat length is "<<mMinRepeatLength<<" >= "<< actual_repeat_length <<" <= "<<mMaxRepeatLength, 8);
#endif
                
                if (qcFoundRepeats(read1) && qcFoundRepeats(read2))
                {
                    
#ifdef DEBUG
                    logInfo("Passed all tests!", 8);
                    read1.inspect();
                    logInfo("-------------------", 8)
#endif
                    
                    return true;
                }
            }
#ifdef DEBUG
            else
            {
                logInfo("\tFailed test 2. Repeat length: "<<actual_repeat_length/* << " : " << match_found*/, 8);
            }
#endif
            j = (*(read1.repeatBegin())).second;
        }
        read1.clearRepeatPositions();
        read2.clearRepeatPositions();
    }
    return false;
}

gzFile crass::Search::getFileHandle(const char *inputFile) {
    gzFile fp;
    if ( strcmp(inputFile, "-") == 0 )
    {
        fp = gzdopen(fileno(stdin), "r");
    }
    else
    {
        fp = gzopen(inputFile, "r");
    }

    if ( (fp == NULL) && (strcmp(inputFile, "-") != 0) )
    {
        std::cerr<<"[ERROR] Could not open FASTQ "<<inputFile<<" for reading."<<std::endl;
        exit(1);
    }

    if ( (fp == NULL) && (strcmp(inputFile, "-") == 0) )
    {
        std::cerr<<"[ERROR] Could not open stdin for reading."<<std::endl;
        exit(1);
    }
    return fp;
}

int crass::Search::scanRight(crass::RawRead& read,
                             std::string& pattern)
{
#ifdef DEBUG
    logInfo("Scanning Right for more repeats:", 9);
#endif
    
    unsigned int pattern_length = static_cast<int>(pattern.length());
    
    // final start index
    auto it = read.repeatEnd();
    int last_repeat_index = (*(it - 1)).first;
    
    //second to final start index
    int second_last_repeat_index = (*(it - 2)).first;
    
    int repeat_spacing = last_repeat_index - second_last_repeat_index;
    
#ifdef DEBUG
    logInfo(read.numberOfRepeats()<<" : "<<pattern_length<<" : "<<last_repeat_index<<" : "<<second_last_repeat_index<<" : "<<repeat_spacing, 9);
#endif
    
    int candidate_repeat_index, position;
    
    int begin_search, end_search;
    
    int read_length = static_cast<int>(read.length());
    bool more_to_search = true;
    while (more_to_search)
    {
        candidate_repeat_index = last_repeat_index + repeat_spacing;
        begin_search = candidate_repeat_index - mScanRange;
        end_search = candidate_repeat_index + pattern_length + mScanRange;
#ifdef DEBUG
        logInfo(candidate_repeat_index<<" : "<<begin_search<<" : "<<end_search, 9);
#endif
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right
        int scanRightMinBegin = last_repeat_index + pattern_length + mMinSpacerLength;
        
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
        
        std::string text = read.seq().substr(begin_search, (end_search - begin_search));
        
#ifdef DEBUG
        logInfo(pattern<<" : "<<text, 9);
#endif
        position = PatternMatcher::bmpSearch(text, pattern);
        
        
        if (position >= 0)
        {
            logInfo(begin_search<<" : "<<position<<" : "<<pattern_length, 10);
            read.push_back(begin_search + position, begin_search + position + pattern_length);
            second_last_repeat_index = last_repeat_index;
            last_repeat_index = begin_search + position;
            repeat_spacing = last_repeat_index - second_last_repeat_index;
            if (repeat_spacing < (mMinSpacerLength + pattern_length))
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


int crass::Search::extendPreRepeat(crass::RawRead& read)
{
    
    int num_repeats = read.numberOfRepeats();
    int cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * num_repeats);
    
    // make sure that we don't go below 2 or everything will match
    if (2 > cut_off)
    {
        cut_off = 2;
    }
#ifdef DEBUG
    read.inspect();
    logInfo("cutoff: "<<cut_off, 9);
#endif
    RepeatArray::RepeatIterator it = read.repeatBegin();
    int first_repeat_start_index = (*it).first;
    
    it = read.repeatEnd() - 1;
    int last_repeat_start_index = (*it).first;
    

    int shortest_repeat_spacing = (*(read.repeatAt(1))).first - (*(read.repeatAt(0))).first;

    RepeatArray::RepeatIterator it2;
    for (it = read.repeatBegin(), it2 = it + 1; it2 != read.repeatEnd(); ++it, ++it2) {
        int curr_repeat_spacing = (*it2).first - (*it).first;
        
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }

    int right_extension_length = 0;
    int max_right_extension_length = shortest_repeat_spacing - mMinSpacerLength;

#ifdef DEBUG
    logInfo("max right extension length: "<<max_right_extension_length, 9);
    logInfo("shortest repeat spacing: "<<shortest_repeat_spacing, 9);

#endif
    RepeatArray::RepeatIterator end_iter = read.repeatEnd();
    RepeatArray::RepeatIterator begin_iter = read.repeatBegin();
    
    
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;

    while (max_right_extension_length > right_extension_length )
    {
        if((last_repeat_start_index + mSearchWindowLength + right_extension_length) >= static_cast<int>(read.length()))
        {
            --end_iter;
        }
        for (RepeatArray::RepeatIterator it = read.repeatBegin(); it != end_iter; ++it )
        {
            switch( read[(*it).second + right_extension_length])
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
        logInfo("R: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << right_extension_length, 9);
#endif
        if ( (char_count_A >= cut_off) || (char_count_C >= cut_off) || (char_count_G >= cut_off) || (char_count_T >= cut_off) )
        {
            right_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }
    
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    
    int left_extension_length = 0;
    int max_left_extension_length = shortest_repeat_spacing - mMinSpacerLength - right_extension_length;
    if(max_left_extension_length < 0) {
        max_left_extension_length = 0;
    }
#ifdef DEBUG
    logInfo("max left extension length: "<<max_left_extension_length, 9);
#endif
    
    begin_iter = read.repeatBegin();
    end_iter = read.repeatEnd();
    while (left_extension_length <= max_left_extension_length)
    {
        if(first_repeat_start_index - left_extension_length < 0) {
            logInfo("first repeat can't be extended anymore, dorpping", 10);
            ++begin_iter;
        }
        for (RepeatArray::RepeatIterator it = begin_iter; it != read.repeatEnd(); ++it )
        {
            switch(read[(*it).first - left_extension_length - 1])
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
        logInfo("L:" << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : "  << " : " << left_extension_length, 9);
#endif
        
        if ( (char_count_A >= cut_off) || (char_count_C >= cut_off) || (char_count_G >= cut_off) || (char_count_T >= cut_off) )
        {
            left_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }

    for(auto it = read.repeatBegin(); it != read.repeatEnd(); ++it) {
        if((*it).first - left_extension_length < 0){
            (*it).first = 0;
        } else {
            (*it).first -= left_extension_length;
        }
        
        if((*it).second + right_extension_length >= static_cast<int>(read.length())) {
            (*it).second = static_cast<int>(read.length());
        } else {
            (*it).second += right_extension_length;
        }
    }
    return read.getFirstNonPartialRepeatLength();
}

int crass::Search::extendPreRepeatPair(crass::RawRead& templateRead, crass::RawRead& read2)
{
    
    int num_repeats = templateRead.numberOfRepeats() + read2.numberOfRepeats();

    
    RepeatArray::RepeatIterator it = templateRead.repeatBegin();
    RepeatArray::RepeatIterator it2 = read2.repeatBegin();

    int template_read_first_repeat_start_index = (*it).first;
    int read2_first_repeat_start_index = (*it2).first;
    
    it = templateRead.repeatEnd() - 1;
    it2 = read2.repeatEnd() - 1;
    int template_read_last_repeat_start_index = (*it).first;
    int read2_last_repeat_start_index = (*it2).first;

    
    int shortest_repeat_spacing = (*(templateRead.repeatAt(1))).first - (*(templateRead.repeatAt(0))).first;

    for (it = templateRead.repeatBegin(), it2 = it + 1; it2 != templateRead.repeatEnd(); ++it, ++it2) {
        int curr_repeat_spacing = (*it2).first - (*it).first;
        
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }
    if (read2.numberOfRepeats() > 1) {
        for (it = read2.repeatBegin(), it2 = it + 1; it2 != read2.repeatEnd(); ++it, ++it2) {
            int curr_repeat_spacing = (*it2).first - (*it).first;
            
            if (curr_repeat_spacing < shortest_repeat_spacing)
            {
                shortest_repeat_spacing = curr_repeat_spacing;
            }
        }
    }



    int right_extension_length = 0;
    int max_right_extension_length = shortest_repeat_spacing - mMinSpacerLength;
    int left_extension_length = 0;
    int max_left_extension_length = shortest_repeat_spacing - mMinSpacerLength - right_extension_length;
    if(max_left_extension_length < 0) {
        max_left_extension_length = 0;
    }

    
#ifdef DEBUG
    logInfo("shortest repeat spacing: "<<shortest_repeat_spacing, 9);
    logInfo("max right extension: "<<max_right_extension_length, 9);
    logInfo("max left extension: "<<max_left_extension_length, 9);
#endif
    RepeatArray::RepeatIterator template_read_end_iter = templateRead.repeatEnd();
    RepeatArray::RepeatIterator read2_end_iter = read2.repeatEnd();
    RepeatArray::RepeatIterator template_read_begin_iter = templateRead.repeatBegin();
    RepeatArray::RepeatIterator read2_begin_iter = read2.repeatBegin();

    
    int cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * num_repeats);
    
    if (2 > cut_off)
    {
        cut_off = 2;
    }
#ifdef DEBUG
    logInfo("cuoff: "<<cut_off, 9);
#endif
    
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;

    while (max_right_extension_length > right_extension_length)
    {
        for ( it = templateRead.repeatBegin(); it != template_read_end_iter; ++it )
        {

            if(((*it).second + right_extension_length) >= static_cast<int>(templateRead.length()))
            {
                it = template_read_end_iter;
            }
            else
            {
                switch( templateRead[(*it).second + right_extension_length])
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
        
        for ( it = read2.repeatBegin(); it != read2_end_iter; ++it )
        {
            
            if(((*it).second + right_extension_length) >= static_cast<int>(templateRead.length()))
            {
                it = read2_end_iter;
            }
            else
            {
                switch( read2[(*it).second + right_extension_length])
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
        logInfo("R: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << right_extension_length, 9);
#endif
        if ( (char_count_A >= cut_off) || (char_count_C >= cut_off) || (char_count_G >= cut_off) || (char_count_T >= cut_off) )
        {
            right_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }
    
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    

    while (left_extension_length <= max_left_extension_length)
    {
        for ( it = template_read_begin_iter; it != templateRead.repeatEnd(); ++it )
        {
            switch(templateRead[(*it).first - left_extension_length - 1])
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
        
        for ( it = read2_begin_iter; it != read2.repeatEnd(); ++it )
        {
            switch(read2[(*it).first - left_extension_length - 1])
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
        logInfo("L: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : "  << left_extension_length, 9);
#endif
        if ( (char_count_A >= cut_off) || (char_count_C >= cut_off) || (char_count_G >= cut_off) || (char_count_T >= cut_off) )
        {
            left_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            break;
        }
    }
    
    for( it = templateRead.repeatBegin(); it != templateRead.repeatEnd(); ++it) {
        if((*it).first - left_extension_length < 0){
            (*it).first = 0;
        } else {
            (*it).first -= left_extension_length;
        }
        
        if((*it).second + right_extension_length >= static_cast<int>(templateRead.length())) {
            (*it).second = static_cast<int>(templateRead.length());
        } else {
            (*it).second += right_extension_length;
        }
    }
    
    for( it = read2.repeatBegin(); it != read2.repeatEnd(); ++it) {
        if((*it).first - left_extension_length < 0){
            (*it).first = 0;
        } else {
            (*it).first -= left_extension_length;
        }
        
        if((*it).second + right_extension_length >= static_cast<int>(read2.length())) {
            (*it).second = static_cast<int>(read2.length());
        } else {
            (*it).second += right_extension_length;
        }
    }
    return templateRead.getFirstNonPartialRepeatLength();
}

bool crass::Search::qcFoundRepeats(crass::RawRead& read)
{

    std::string repeat = *(read.getFirstNonPartialRepeat());
    
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
    if (read.numberOfRepeats() == 1) {
        return true;
    }
    // test for a long or short read
    int single_compare_index = 0;
    bool is_short = (2 > read.numberOfSpacers());
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
        
        // now go through the spacers and check for similarities
        std::vector<std::string> spacer_vec;
        
        //get all the spacer strings
        for (auto it = read.spacerStringBegin(); it != read.spacerStringEnd(); ++it) {
            spacer_vec.push_back(*it);
        }
        
        std::vector<std::string>::iterator spacer_iter = spacer_vec.begin();
        // comparing pairs below so minus 1 from end
        std::vector<std::string>::iterator spacer_last = spacer_vec.end() - 1;
        
        
        for ( ;spacer_iter != spacer_last; ++spacer_iter)
        {
            
            num_compared++;
            ave_repeat_to_spacer_difference += PatternMatcher::getStringSimilarity(repeat, *spacer_iter);

            float ss_diff = 0;
            ss_diff += PatternMatcher::getStringSimilarity(*spacer_iter, *(spacer_iter + 1));

            ave_spacer_to_spacer_difference += ss_diff;
            
            ave_spacer_to_spacer_len_difference += (static_cast<float>(spacer_iter->size()) - static_cast<float>((spacer_iter + 1)->size()));
            ave_repeat_to_spacer_len_difference +=  (static_cast<float>(repeat.size()) - static_cast<float>(spacer_iter->size()));
        }
        
        // now look for max and min lengths!
        spacer_iter = spacer_vec.begin();
        spacer_last = spacer_vec.end();
        while (spacer_iter != spacer_last)
        {
            if(static_cast<int>(spacer_iter->length()) < min_spacer_length)
            {
                min_spacer_length = static_cast<int>(spacer_iter->length());
            }
            if(static_cast<int>(spacer_iter->length()) > max_spacer_length)
            {
                max_spacer_length = static_cast<int>(spacer_iter->length());
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
            ave_spacer_to_spacer_difference /= static_cast<float>(num_compared);
            ave_repeat_to_spacer_difference /= static_cast<float>(num_compared);
            ave_spacer_to_spacer_len_difference = abs(ave_spacer_to_spacer_len_difference /= static_cast<float>(num_compared));
            ave_repeat_to_spacer_len_difference = abs(ave_repeat_to_spacer_len_difference /= static_cast<float>(num_compared));
            
            /*
             * MAX AND MIN SPACER LENGTHS
             */
            if (min_spacer_length < mMinSpacerLength)
            {
#ifdef DEBUG
                logInfo("\tFailed test 4a. Min spacer length out of range: "<<min_spacer_length<<" < "<<mMinSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4a. Min spacer length within range: "<<min_spacer_length<<" > "<<mMinSpacerLength, 8);
#endif
            if (max_spacer_length > mMaxSpacerLength)
            {
#ifdef DEBUG
                logInfo("\tFailed test 4b. Max spacer length out of range: "<<max_spacer_length<<" > "<<mMaxSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4b. Max spacer length within range: "<<max_spacer_length<<" < "<<mMaxSpacerLength, 8);
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
        std::string spacer = *(read.spacerStringAt(single_compare_index));
        float similarity = PatternMatcher::getStringSimilarity(repeat, spacer);
        if (similarity > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY)
        {
            /*
             * MAX AND MIN SPACER LENGTHS
             */
            if (static_cast<int>(spacer.length()) < mMinSpacerLength)
            {
#ifdef DEBUG
                logInfo("\tFailed test 4a. Min spacer length out of range: "<<spacer.length()<<" < "<<mMinSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4a. Min spacer length within range: "<<spacer.length()<<" > "<<mMinSpacerLength, 8);
#endif
            if (static_cast<int>(spacer.length()) > mMaxSpacerLength)
            {
#ifdef DEBUG
                logInfo("\tFailed test 4b. Max spacer length out of range: "<<spacer.length()<<" > "<<mMaxSpacerLength, 8);
#endif
                return false;
            }
#ifdef DEBUG
            logInfo("\tPassed test 4b. Max spacer length within range: "<<spacer.length()<<" < "<<mMaxSpacerLength, 8);
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
        if (abs(static_cast<int>(spacer.length()) - static_cast<int>(repeat.length())) > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF)
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

bool crass::Search::isRepeatLowComplexity(std::string& repeat)
{
    int c_count = 0;
    int g_count = 0;
    int a_count = 0;
    int t_count = 0;
    int n_count = 0;
    
    int curr_repeat_length = static_cast<int>(repeat.length());
    
    int cut_off = static_cast<int>(curr_repeat_length * CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD);
    
    std::string::iterator dr_iter;
    for (dr_iter = repeat.begin(); dr_iter != repeat.end();++dr_iter)
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
    }
    if (a_count > cut_off) return true;
    else if (t_count > cut_off) return true;
    else if (g_count > cut_off) return true;
    else if (c_count > cut_off) return true;
    else if (n_count > cut_off) return true;
    return false;
}


bool crass::Search::doesRepeatHaveHighlyAbundantKmers(std::string& directRepeat) {
    float max;
    return doesRepeatHaveHighlyAbundantKmers(directRepeat, max);
}

bool crass::Search::doesRepeatHaveHighlyAbundantKmers(std::string& directRepeat, float& maxFrequency)
{
    // cut kmers from the direct repeat to test whether
    // a particular kmer is vastly over represented
    std::map<std::string, int> kmer_counter;
    size_t kmer_length = 3;
    size_t max_index = (directRepeat.length() - kmer_length);
    std::string kmer;
    int total_count = 0;
        for (size_t i = 0; i < max_index; i++) {
            kmer = directRepeat.substr(i, kmer_length);
            if (kmer_counter.find(kmer) == kmer_counter.end()) {
                kmer_counter[kmer] = 1;
            } else {
                kmer_counter[kmer]++;
            }
            total_count++;
        }
    
    std::map<std::string, int>::iterator iter;
    int max_count = 0;
    for (iter = kmer_counter.begin(); iter != kmer_counter.end(); ++iter) {
        if (iter->second > max_count) {
            max_count = iter->second;
        }
    }
    maxFrequency = static_cast<float>(max_count)/static_cast<float>(total_count);
    if (maxFrequency > CRASS_DEF_KMER_MAX_ABUNDANCE_CUTOFF) {
        return true;
    } else {
        return false;
    }
}

#if crass_Search_main
int main(int argc, char *argv[]) {
    crass::Search s = crass::Search();
    s.minSeedCount(2);
    //                                                                                                                                                       1                                                                                                   2                                                                                                   3
    //                                                             1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4
    //                                                   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    //crass::RawRead read100("Accumulibacter_read100", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATG", "", crass::RepeatArray());
    //crass::RawRead read350("Accumulibacter_CRISPR", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATGGCCTGGCGCATGGGCTTGGCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCCCATCGCTCCAGGCGACACGCCGCAAGCGGCGGTTTCCCCCGCGTCAGCGGGGATAGGCCCACGCGGCCATCGAATTCGGTGGTCTTGCCGCCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTGTGCTGGCGCGTGAAGTCGCAGCCCATTGGCGTTTCCCCCGCGTCAGCGGGGATAGGCCGCAAAGCCACAATCTTT", "",crass::RepeatArray());
    //Read350                                            RRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSS
    //Read100                                            RRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSS
    
    //crass::RawRead read100_mate("Accumulibacter_read100_mate", "", "CGTCAGCGGGGATAGGCCTCAGCCCCTACGAGCGCACCGACGACCCGATTGGTTTCCCCCGCGTCAGCGGGGATAGGCCCTTTGGGCGCAAGGATGACGC", "", crass::RepeatArray());
    
    //s.longReadSearch(read350);
    //s.longReadSearch(read100);
    //s.longReadPairSearch(read100, read100_mate);
    //read350.inspect();
    //read100.inspect();
    //read100_mate.inspect();
    
    //crass::RawRead read1_not("", "", "AGGTAGTGTCGTCGTTGTCTGGCGCGATCTGAATGAAACCGGGGGTGTTGGGATCAAGACGCTGGGTGAGCCGTCGAAGGAGCTTGTCGAGGTTGACGGCCTCTGGCTGGTGAAACGGATGGTATAAAGTGCTCTTTAACAAAATGAATCGATGTTGAGGTTTTCTCTTGTAGTTCATAGGGCTACCTACAAGAGTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATGGCCTGG", "", crass::RepeatArray());
    //crass::RawRead read2_yes("", "", "ACCCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCCGGCGTGTTCGAATTCAATGGTCTTGCGCTCGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCCTGCCGGCTGGCCTTGCGCATCGGGCTCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGGGCGGCCCATACCGGTGACGCCCTTGGTAGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCCTTCCCGGCGCATGTCGTCGGCACCTGACTAGGTTTCCCCCGCGTCAGCGGGGATAGGCCCTTTTCGATCACTTCCTGGCCT", "", crass::RepeatArray());
    
    //s.longReadPairSearch(read1_not, read2_yes);
    //read1_not.inspect();
    //read2_yes.inspect();
    //std::cout<<"*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*"<<std::endl;
    //crass::RawRead read1_yes("", "", "GGCAGTTTCCCCCGCGTCAGCGGGGATAGGCCCCATCGACGACGGCCAGCGCCTGGTCGAACAGTGTTTCCCCCGCGTCAGCGGGGATAGGCCCGGTGGGCGAGATTCGCGACTCGCGGCCGTTGCGTTTCCCCCGCGTCAGCGGGGATAGGCCCCGCAGCACGGGCAGGGCTCGCAGCGCTGCGAGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGTCGTCCGGCTTGTGGGATTCGGCAGTATGGGGTTTCCCCCGCGTCAGCGGGGATAGGCCGCCAAGGGTCGAGGGCTTGGCGT", "", crass::RepeatArray());
    //crass::RawRead read2_not("", "", "GGCGACAGTCTCGCCCCCGGCTGACCTACAGTTTCCCCCGCGTCAGCGGGGATAGGCCCCGGGCGTGCCCGAAGTCGTCGAGAGATGCTTAGGCGCTGCGGGGCCTGCCGCGTTCGCTGTAGTAGATCGCGGCCTCGTCGAGTTCGGCCAGCGCCGGCGGCAGGATGACGATTTTCATGTCGGCCGGTGGCTGGCGATTCTGGCGCGGGTCTGCGCAATGGCTTCGTCGATATCAACGGTGCCAATCTCACCGCGGTCGAAAGCGTCGGCTCTGCGCTCGAGTTCTTCCTGCCACGCGGCCAGG", "", crass::RepeatArray());
    //s.longReadSearch(read1_yes);
    //s.longReadPairSearch(read1_yes, read2_not);
    //s.longReadPairSearch(read2_not, read1_yes);
    //read1_yes.inspect();
    //read2_not.inspect();
    
    //crass::RawRead read1_bothyes("", "", "AAAGGGTTTGGCGACGCAAGTTTCCCCCGCGTCAGCGGGGATAGGCCCCTCGAGCTGGACCGCAATCTGGAGTCTGCCCTGTTTCCCCCGCGTCAGCGGGGATAGGCCCGACGCCATCGAGGCGATCGAGGCGCAGATCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGTCCGAGCATCGCGGCGTATCAAAGCGCGGGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCTCCACAGCGTAGAAATGCCG", "", crass::RepeatArray());
    //crass::RawRead read2_bothyes("", "", /*crass::reverseComplement(*/"TGTTTCCCCCGCGTCAGCGGGGATAGGCCCCCGCACCACGACGACACCAGCTTACTACGGCCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTCGCCCATATCGGTATCTGTCCAGAGCGCAGCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTGCTGCGGGATCATCGCGTCGGCACGCCCGGCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTGTACGGTCGGGGGTGCTGGATATATATCGAGGTTT"/*)*/, "", crass::RepeatArray());
    //s.longReadSearch(read1_bothyes);
    //s.longReadSearch(read2_bothyes);
    //s.longReadPairSearch(read1_bothyes, read2_bothyes);
    //read1_bothyes.inspect();
    //read2_bothyes.inspect();
    
    //crass::RawRead trap_read("trap", "", "CCCGCACGCATTAAGGCGGTTCATCCCCGCGCCTGCGGGGAACGCAAATTGATAGCGCCAAGCGCCACAGCTTGCGCCGGTTCATCCCCGCGCCTGCGGGG", "", crass::RepeatArray());
    //s.longReadSearch(trap_read);
    //trap_read.inspect();
    if(argc >1) {
        s.searchFileSerial(argv[1]);
    }
    
    return 0;
}
#endif

