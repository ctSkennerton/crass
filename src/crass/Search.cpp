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

// local includes
#include "Search.h"
#include "PatternMatcher.h"



int longReadSearch(ReadHolder& tmpHolder,
                   const options& opts,
                   ReadMap * mReads,
                   StringCheck * mStringCheck,
                   lookupTable& patternsHash,
                   lookupTable& readsFound)
{
    //-----
    // Code lifted from CRT, ported by Connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    
    //bool match_found = false;
    
    std::string read = tmpHolder.getSeq();
    
    // get the length of this sequence
    unsigned int seq_length = static_cast<unsigned int>(read.length());
    
    
    //the number of bases that can be skipped while we still guarantee that the entire search
    //window will at some point in its iteration thru the sequence will not miss a any repeat
    unsigned int skips = mMinRepeatLength - (2 * mSearchWindowLength - 1);
    if (skips < 1)
    {
        skips = 1;
    }
    
    int searchEnd = seq_length - opts.highDRsize - opts.highSpacerSize - mSearchWindowLength - 1;
    
    if (searchEnd < 0)
    {
        logError("Read: "<<tmpHolder.getHeader()<<" length is less than "<<opts.highDRsize + opts.highSpacerSize + mSearchWindowLength + 1<<"bp");
        //delete tmpHolder;
        return 1;
    }
    
    for (unsigned int j = 0; j <= static_cast<unsigned int>(searchEnd); j = j + skips)
    {
        
        unsigned int beginSearch = j + mMinRepeatLength + opts.lowSpacerSize;
        unsigned int endSearch = j + opts.highDRsize + opts.highSpacerSize + mSearchWindowLength;
        
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
            throw crispr::substring_exception(e.what(),
                                              text.c_str(),
                                              beginSearch,
                                              (endSearch - beginSearch),
                                              __FILE__,
                                              __LINE__,
                                              __PRETTY_FUNCTION__);
        }
        try {
            pattern = read.substr(j, mSearchWindowLength);
        } catch (std::exception& e) {
            throw crispr::substring_exception(e.what(),
                                              read.c_str(),
                                              j,
                                              mSearchWindowLength,
                                              __FILE__,
                                              __LINE__,
                                              __PRETTY_FUNCTION__);
        }
        
        //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
        int pattern_in_text_index = -1;
        pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);
        
        if (pattern_in_text_index >= 0)
        {
            tmpHolder.startStopsAdd(j,  j + mSearchWindowLength);
            unsigned int found_pattern_start_index = beginSearch + static_cast<unsigned int>(pattern_in_text_index);
            
            tmpHolder.startStopsAdd(found_pattern_start_index, found_pattern_start_index + mSearchWindowLength);
            scanRight(tmpHolder, pattern, opts.lowSpacerSize, 24);
        }
        
        if ( (tmpHolder.numRepeats() >= opts.minNumRepeats) ) //tmp_holder->numRepeats is half the size of the StartStopList
        {
#ifdef DEBUG
            logInfo(tmpHolder.getHeader(), 8);
            logInfo("\tPassed test 1. At least "<<opts.minNumRepeats<< " ("<<tmpHolder.numRepeats()<<") repeated kmers found", 8);
#endif
            
            unsigned int actual_repeat_length = extendPreRepeat(tmpHolder, mSearchWindowLength, opts.lowSpacerSize);
            
            if ( (actual_repeat_length >= mMinRepeatLength) && (actual_repeat_length <= opts.highDRsize) )
            {
#ifdef DEBUG
                
                logInfo("\tPassed test 2. The repeat length is "<<mMinRepeatLength<<" >= "<< actual_repeat_length <<" <= "<<opts.highDRsize, 8);
#endif
                // Declare a tmp string here to hold the encoded DR if
                // removeHolopolymers is in affect.  Later if the read passes all
                // the tests then add in this encoded string since that is the
                // version that the singleton finder should be looking for
                std::string encoded_repeat;
                if(opts.removeHomopolymers) {
                    encoded_repeat = tmpHolder.repeatStringAt(0);
                    tmpHolder.decode();
                    
                }
                
                // drop partials
                tmpHolder.dropPartials();
                if (qcFoundRepeats(tmpHolder, opts.lowSpacerSize, opts.highSpacerSize))
                {
#ifdef DEBUG
                    logInfo("Passed all tests!", 8);
                    logInfo("Potential CRISPR containing read found: "<<tmpHolder.getHeader(), 7);
                    logInfo(tmpHolder.getSeq(), 9);
                    logInfo("-------------------", 7)
#endif
                    
                    //ReadHolder * candidate_read = new ReadHolder();
                    //*candidate_read = *tmp_holder;
                    //addReadHolder(mReads, mStringCheck, candidate_read);
                    addReadHolder(mReads, mStringCheck, tmpHolder);
                    //match_found = true;
					if(opts.removeHomopolymers) {
                        patternsHash[encoded_repeat] = true;
                    } else {
                        patternsHash[tmpHolder.repeatStringAt(0)] = true;
                    }
                    readsFound[tmpHolder.getHeader()] = true;
                    break;
                }
            }
#ifdef DEBUG
            else
            {
                logInfo("\tFailed test 2. Repeat length: "<<tmpHolder.getRepeatLength()/* << " : " << match_found*/, 8);
            }
#endif
            j = tmpHolder.back() - 1;
        }
        tmpHolder.clearStartStops();
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

#if crass_Search_main
int main() {
    crass::Search s = crass::Search();
    //                                                                                                                                                       1                                                                                                   2                                                                                                   3
    //                                                             1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4
    //                                                   01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    crass::RawRead read100("Accumulibacter_CRISPR", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATG", "", crass::RepeatArray());
    crass::RawRead read350("Accumulibacter_CRISPR", "", "GTTTCCCCCGCGTCAGCGGGGATAGGCCCCACGCCGTGACGGAAGGGCTGATCACGAAATGGTTTCCCCCGCGTCAGCGGGGATAGGCCCGCGCGGCATGGCCTGGCGCATGGGCTTGGCGGGTTTCCCCCGCGTCAGCGGGGATAGGCCCCCATCGCTCCAGGCGACACGCCGCAAGCGGCGGTTTCCCCCGCGTCAGCGGGGATAGGCCCACGCGGCCATCGAATTCGGTGGTCTTGCCGCCGTTTCCCCCGCGTCAGCGGGGATAGGCCCTGTGCTGGCGCGTGAAGTCGCAGCCCATTGGCGTTTCCCCCGCGTCAGCGGGGATAGGCCGCAAAGCCACAATCTTT", "",crass::RepeatArray());
    //                                                   RRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    
    read350.push_back(0,6);
    read350.push_back(61,70);
    std::string str = "GTTTCCCCC";
    s.scanRight(read350, str);
    int new_length = s.extendPreRepeat(read350);
    std::cout<<"Extended Length: "<< new_length<<std::endl;
    read350.inspect();
    return 0;
}
#endif

