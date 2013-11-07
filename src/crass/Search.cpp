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

int crass::Search::shortReadSearch(crass::RawRead&  read)
{
    
    //bool match_found = false;
    
    
    int seq_length = static_cast<int>(read.seq().length());
    int search_end = seq_length - mMinRepeatLength - 1;
    int final_index = seq_length - 1;
    
    
    // boyer-moore search
    for (int first_start = 0; first_start < search_end; first_start++)
    {
        int search_begin = first_start + mMinRepeatLength + mMinSpacerLength;
        
        if (search_begin >= search_end ) break;
        
        // do the search
        int second_start = -1;
        
        //try {
            second_start = PatternMatcher::bmpSearch(read.seq().substr(search_begin),
                                                     read.seq().substr(first_start, mMinRepeatLength) );
        /*} catch (std::exception& e) {
            
			std::stringstream ss;
            ss<<read<<"\n";
            ss<<search_begin<<" : "<<first_start<<" : "<<opts.lowDRsize;
            throw crispr::runtime_exception( __FILE__,
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            ss);
        }*/
        
        // check to see if we found something
        if (second_start > -1)
        {
            // bingo!
            second_start += search_begin;
            int second_end = second_start + mMinRepeatLength;
            int first_end = first_start + mMinRepeatLength;
            
            int next_index = second_end + 1;
            
            // make sure that the kmer match is not already at the end of the read before incrementing
            // increment so we are looking at the next base after the match
            if ( next_index <= final_index)
            {
                // read through the subsuquent bases untill they don't match
                int extenstion_length = 0;
                while (read[first_end + extenstion_length] == read[second_end + extenstion_length])
                {
#ifdef DEBUG
                    logInfo(read[first_end + extenstion_length]<<" == "<<read[second_end + extenstion_length],20);
#endif
                    extenstion_length++;
                    next_index++;
                    if (next_index > final_index) break;
                    
                }
#ifdef DEBUG
                logInfo("A: FS: " << first_start<< " FE: "<<first_end<<" SS: "<< second_start<< " SE: "<<second_end<<" EX: "<<extenstion_length, 8);
#endif
                read.push_back(first_start, first_end + extenstion_length);
                read.push_back(second_start, second_end + extenstion_length);
            }
            else
            {
#ifdef DEBUG
                logInfo("B: FS: " << first_start<< " FE: "<<first_end<<" SS: "<< second_start<< " SE: "<<second_end, 8);
#endif
                read.push_back(first_start, first_end);
                read.push_back(second_start, second_end);
            }
			
            // the low side will always be true since we search for the lowDRSize
            if ( read.getFirstNonPartialRepeatLength() <= mMaxRepeatLength )
            {
                int i = static_cast<int>((*read.spacerStringAt(0)).length());
                if (( i >= mMinSpacerLength) && (i <= mMaxSpacerLength))
                {
                    if (qcFoundRepeats(read))
                    {
                        
                        //match_found = true;
#ifdef DEBUG
                        logInfo("Potential CRISPR containing read found: ", 8);
                        read.inspect();
                        logInfo("-------------------", 8)
#endif
                        //readsFound[tmpHolder.getHeader()] = true;
                        //addReadHolder(mReads, mStringCheck, tmpHolder);
                        break;
                    }
                }
#ifdef DEBUG
                else
                {
                    logInfo("\tFailed test 2. The spacer length is not between "<<mMinSpacerLength<<" and "<<mMaxSpacerLength<<": "<< i, 8);
                }
#endif
            }
#ifdef DEBUG
            else
            {
                logInfo("\tFailed test 1. The repeat length is larger than "<<mMaxRepeatLength<<": " << read.getFirstNonPartialRepeatLength(), 8);
            }
#endif
            first_start = (*(read.repeatEnd() - 1)).second;
        }
    }
    return 0;
}


int crass::Search::longReadSearch(crass::RawRead& read)
{
    //-----
    // Code lifted from CRT, ported by Connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    
    //bool match_found = false;
        
    // get the length of this sequence
    int seq_length = static_cast<int>(read.seq().length());
    
    
    //the number of bases that can be skipped while we still guarantee that the entire search
    //window will at some point in its iteration thru the sequence will not miss a any repeat
    unsigned int skips = mMinRepeatLength - (2 * mSearchWindowLength - 1);
    if (skips < 1)
    {
        skips = 1;
    }
    
    int searchEnd = seq_length - mMaxRepeatLength - mMaxSpacerLength - mSearchWindowLength;
    
    if (searchEnd < 0)
    {
        //logError("Read:"<<read.identifier()<<" length is less than " << mMaxRepeatLength + mMaxSpacerLength + mSearchWindowLength<<"bp");
        //delete tmpHolder;
        return 1;
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
        //try {
            text = read.seq().substr(beginSearch, (endSearch - beginSearch));
            
        /*} catch (std::exception& e) {
            throw crispr::substring_exception(e.what(),
                                              text.c_str(),
                                              beginSearch,
                                              (endSearch - beginSearch),
                                              __FILE__,
                                              __LINE__,
                                              __PRETTY_FUNCTION__);
        }
        try {*/
            pattern = read.seq().substr(j, mSearchWindowLength);
        /*} catch (std::exception& e) {
            throw crispr::substring_exception(e.what(),
                                              read.c_str(),
                                              j,
                                              mSearchWindowLength,
                                              __FILE__,
                                              __LINE__,
                                              __PRETTY_FUNCTION__);
        }*/
        
        //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
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
                
                // drop partials
                //tmpHolder.dropPartials();
                if (qcFoundRepeats(read))
                {
#ifdef DEBUG
                    logInfo("Passed all tests!", 8);
                    read.inspect();
                    logInfo("-------------------", 8)
#endif
                    
                    //ReadHolder * candidate_read = new ReadHolder();
                    //*candidate_read = *tmp_holder;
                    //addReadHolder(mReads, mStringCheck, candidate_read);
                    //addReadHolder(mReads, mStringCheck, tmpHolder);
                    //match_found = true;
//					if(opts.removeHomopolymers) {
//                        patternsHash[encoded_repeat] = true;
//                    } else {
//                        patternsHash[tmpHolder.repeatStringAt(0)] = true;
//                    }
//                    readsFound[tmpHolder.getHeader()] = true;
                    break;
                }
            }
#ifdef DEBUG
            else
            {
                logInfo("\tFailed test 2. Repeat length: "<<actual_repeat_length/* << " : " << match_found*/, 8);
            }
#endif
            j = (*(read.repeatEnd() - 1)).second;
        }
        read.clearRepeatPositions();
    }
    return 0;
}

int crass::Search::scanRight(crass::RawRead& read,
                             std::string& pattern)
{
#ifdef DEBUG
    logInfo("Scanning Right for more repeats:", 9);
#endif
    //unsigned int start_stops_size = tmp_holder.getStartStopListSize();
    
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
#ifdef DEBUG
    //tmp_holder.logContents(9);
#endif
#ifdef DEBUG
    logInfo("Extending Prerepeat...", 9);
#endif
    
    int num_repeats = read.numberOfRepeats();
    //tmp_holder.setRepeatLength(searchWindowLength);
    int cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * num_repeats);
    
    // make sure that we don't go below 2
    if (2 > cut_off)
    {
        cut_off = 2;
    }
#ifdef DEBUG
    logInfo("cutoff: "<<cut_off, 9);
#endif
    read.inspect();
    RepeatArray::RepeatIterator it = read.repeatBegin();
    // the index in the read of the first DR kmer
    int first_repeat_start_index = (*it).first;
    
    // the index in the read of the last DR kmer
    it = read.repeatEnd() - 1;
    int last_repeat_start_index = (*it).first;
    
    // the length between the first two DR kmers
    int shortest_repeat_spacing = (*(read.repeatAt(1))).first - (*(read.repeatAt(0))).first;
    // loop througth all remaining members of mRepeats
    //unsigned int end_index = tmp_holder.getStartStopListSize();
    RepeatArray::RepeatIterator it2;
    for (it = read.repeatBegin(), it2 = it + 1; it2 != read.repeatEnd(); ++it, ++it2) {
        int curr_repeat_spacing = (*it2).first - (*it).first;
    #ifdef DEBUG
        logInfo("curr repeat spacing : "<<curr_repeat_spacing, 10);
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
    int right_extension_length = 0;
    // don't search too far
    int max_right_extension_length = shortest_repeat_spacing - mMinSpacerLength;
#ifdef DEBUG
    logInfo("max right extension length: "<<max_right_extension_length, 9);
#endif
    // Sometimes we shouldn't use the far right DR. (it may lie too close to the end)
    RepeatArray::RepeatIterator end_iter = read.repeatEnd();
    
    //int DR_index_end = end_index;
    int dist_to_end = read.length() - last_repeat_start_index;
    if(dist_to_end < max_right_extension_length)
    {
#ifdef DEBUG
        logInfo("removing end partial: "<<last_repeat_start_index<<" ("<<dist_to_end<<" < "<<max_right_extension_length<<")", 9);
#endif
        --end_iter;
        //DR_index_end -= 2;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off)
        {
            cut_off = 2;
        }
#ifdef DEBUG
        logInfo("minimum number of repeats needed with same nucleotide for extension: "<<cut_off, 9);
#endif
    }
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    // (from the right side) extend the length of the repeat to the right
    // as long as the last base of all repeats are at least threshold
#ifdef DEBUG
    //int iteration_counter = 0;
#endif
    while (max_right_extension_length > right_extension_length)
    {
        for (auto it = read.repeatBegin(); it != end_iter; ++it )
        {
#ifdef DEBUG
            //logInfo("Iteration "<<iteration_counter++, 10);
            //logInfo(k<<" : "<<tmp_holder.getRepeatAt(k) + tmp_holder.getRepeatLength(), 10);
#endif
            // look at the character just past the end of the last repeat
            // make sure our indicies make some sense!
            if(((*it).second + right_extension_length) >= static_cast<int>(read.length()))
            {
                it = end_iter;
            }
            else
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
        }
        
#ifdef DEBUG
        logInfo("R: " << char_count_A << " : " << char_count_C << " : " << char_count_G << " : " << char_count_T << " : " << max_right_extension_length, 9);
#endif
        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            //tmp_holder.incrementRepeatLength();
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
    int left_extension_length = 0;
    int max_left_extension_length = shortest_repeat_spacing - mMinSpacerLength - right_extension_length;
    if(max_left_extension_length < 0) {
        max_left_extension_length = 0;
    }
#ifdef DEBUG
    logInfo("max left extension length: "<<max_left_extension_length, 9);
#endif
    // and again, we may not wat to use the first DR
    RepeatArray::RepeatIterator begin_iter = read.repeatBegin();
    
    unsigned int DR_index_start = 0;
    if(max_left_extension_length > first_repeat_start_index)
    {
#ifdef DEBUG
        logInfo("removing start partial: "<<max_left_extension_length<<" > "<<first_repeat_start_index, 9);
#endif
        begin_iter++;
        cut_off = static_cast<int>(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
        if (2 > cut_off)
        {
            cut_off = 2;
        }
#ifdef DEBUG
        logInfo("new cutoff: "<<cut_off, 9);
#endif
    }
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (left_extension_length <= max_left_extension_length)
    {
        for (auto it = begin_iter; it != read.repeatEnd(); ++it )
        {
#ifdef DEBUG
            //logInfo(k<<" : "<<tmp_holder.getRepeatAt(k) - left_extension_length - 1, 10);
#endif
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
        
        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            //tmp_holder.incrementRepeatLength();
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
            //try {
            	ave_repeat_to_spacer_difference += PatternMatcher::getStringSimilarity(repeat, *spacer_iter);
            /*} catch (std::exception& e) {
                std::cerr<<"Failed to compare similarity between: "<<repeat<<" : "<<*spacer_iter<<std::endl;
                throw crispr::exception(__FILE__,
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        e.what());
            }*/
            float ss_diff = 0;
            //try {
            	ss_diff += PatternMatcher::getStringSimilarity(*spacer_iter, *(spacer_iter + 1));
			/*} catch (std::exception& e) {
                std::cerr<<"Failed to compare similarity between: "<<*spacer_iter<<" : "<<*(spacer_iter+1)<<std::endl;
                throw crispr::exception(__FILE__,
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        e.what());
			}*/
            ave_spacer_to_spacer_difference += ss_diff;
            
            //MI std::cout << ss_diff << " : " << *spacer_iter << " : " << *(spacer_iter + 1) << std::endl;
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
    //try {
        for (size_t i = 0; i < max_index; i++) {
            kmer = directRepeat.substr(i, kmer_length);
            //std::cout << kmer << ", ";
            if (kmer_counter.find(kmer) == kmer_counter.end()) {
                kmer_counter[kmer] = 1;
            } else {
                kmer_counter[kmer]++;
            }
            total_count++;
        }
   /* } catch (std::exception& e) {
        throw crispr::runtime_exception(__FILE__,
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        e.what());
    }*/
    //std::cout << std::endl;
    
    std::map<std::string, int>::iterator iter;
    //std::cout << directRepeat << std::endl;
    int max_count = 0;
    for (iter = kmer_counter.begin(); iter != kmer_counter.end(); ++iter) {
        if (iter->second > max_count) {
            max_count = iter->second;
        }
        //std::cout << "{" << iter->first << ", " << iter->second << ", " << (float(iter->second)/float(total_count)) << "}, ";
    }
    //std::cout << std::endl;
    maxFrequency = static_cast<float>(max_count)/static_cast<float>(total_count);
    if (maxFrequency > CRASS_DEF_KMER_MAX_ABUNDANCE_CUTOFF) {
        return true;
    } else {
        return false;
    }
}

#if crass_Search_main
int main() {
    crass::Search s = crass::Search();
    //                                                                                                                                                       1                                                                                                   2                                                                                                   3
    //                                                             1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4
    //                                                   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    crass::RawRead read100("Accumulibacter_CRISPR", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATG", "", crass::RepeatArray());
    crass::RawRead read350("Accumulibacter_CRISPR", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATGGCCTGGCGCATGGGCTTGGCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCCCATCGCTCCAGGCGACACGCCGCAAGCGGCGGTTTCCCCCGCGTCAGCGGGGATAGGCCCACGCGGCCATCGAATTCGGTGGTCTTGCCGCCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTGTGCTGGCGCGTGAAGTCGCAGCCCATTGGCGTTTCCCCCGCGTCAGCGGGGATAGGCCGCAAAGCCACAATCTTT", "",crass::RepeatArray());
    //Read350                                            RRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    //Read100                                            RRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSS
    s.longReadSearch(read350);
    read350.inspect();
    
    s.shortReadSearch(read100);
    read100.inspect();
    return 0;
}
#endif

