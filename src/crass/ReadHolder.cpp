// File: ReadHolder.cpp
// Original Author: Michael Imelfort 2011
// Hacked and Extended: Connor Skennerton 2011, 2012
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Implementation of ReadHolder functions
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
#include <sstream>
#include <algorithm>
#include <stdexcept>
// local includes
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "SmithWaterman.h"
#include "LoggerSimp.h"
#include <libcrispr/Exception.h>

unsigned int crass::ReadHolder::getAverageSpacerLength() {
    unsigned int total = 0;
    crispr::RepeatArray<unsigned int>::iterator iter;
    for (iter = RH_StartStops.begin(); iter != RH_StartStops.end() - 1; iter++) {
        total += ((iter+1)->first - 1) - (iter->second + 1);
    }
    return total / RH_StartStops.numberOfSpacers();
}

void crass::ReadHolder::getAllSpacerStrings(std::vector<std::string> &spacers) {
    for (size_t s = 0; s <= RH_StartStops.numberOfSpacers(); s++) {
        spacers.push_back(RH_StartStops.spacerStringAt(s));
    }
}

void crass::ReadHolder::getAllRepeatStrings(std::vector<std::string> &repeats) {
    for (size_t r = 0; r <= RH_StartStops.numberOfRepeats(); r++) {
        repeats.push_back(RH_StartStops.spacerStringAt(r));
    }
}

int crass::ReadHolder::averageRepeatLength() {
    int total = 0;
    crispr::RepeatArray<unsigned int>::iterator iter;
    for (iter = RH_StartStops.begin(); iter != RH_StartStops.end(); iter++) {
        total += iter->second - iter->first + 1;
    }
    return total / RH_StartStops.length();
}

void crass::ReadHolder::startStopsAdd(unsigned int i, unsigned int j) {
#ifdef DEBUG
	if(((int)i < 0) || ((int)j < 0)) {
		std::stringstream ss;
		ss<<"Adding negative to SS list! " << i << " : " << j;
		throw crispr::exception(__FILE__,
		                        __LINE__,
		                        __PRETTY_FUNCTION__,
		                        ss);
	}
	if(i > j) {
		std::stringstream ss;
		ss <<"SS list corrupted! " << i << " : " << j;
		throw crispr::exception(__FILE__,
		                        __LINE__,
		                        __PRETTY_FUNCTION__,
		                        ss);
	}
	if((i > RH_Seq.length()) || (j > RH_Seq.length())) {
		std::stringstream ss;
		ss<<"Too long! " << i << " : " << j;
		throw crispr::exception(__FILE__,
		                        __LINE__,
		                        __PRETTY_FUNCTION__,
		                        ss);
	}
#endif
    if(j >= (unsigned int)getSeqLength())
    {
    	j = (unsigned int)getSeqLength() - 1;
    }
    RH_StartStops.push_back(i, j);
}

void crass::ReadHolder::dropPartials() {
    std::pair<unsigned int, unsigned int> f = RH_StartStops.front();
    if (f.first == 0) {
        RH_StartStops.erase(RH_StartStops.begin());
    }
    f = RH_StartStops.back();
    if (f.second == (unsigned int)getSeqLength() - 1) {
        RH_StartStops.erase(RH_StartStops.end() - 1);
    }
}

void crass::ReadHolder::reverseStartStops() {
    RH_StartStops.reverse();
}

void crass::ReadHolder::updateStartStops(int frontOffset, std::string *DR, const options *opts) {
    //-----
    // Update the start and stops to capture the largest part
    // of the DR
    //
    // Take this opportunity to look for partials at either end of the read
    //
    
    unsigned int DR_length = static_cast<unsigned int>(DR->length());
    
    crispr::RepeatArray<unsigned int>::iterator ss_iter;
    for (ss_iter = RH_StartStops.begin(); ss_iter != RH_StartStops.end(); ss_iter++) {
        unsigned int usable_length = DR_length - 1;
        if (frontOffset >= static_cast<int>(ss_iter->first)) {
            // this will be less than 0
            int amount_below_zero = frontOffset - static_cast<int>(ss_iter->first);
            usable_length = DR_length - amount_below_zero - 1;
            ss_iter->first = 0;
        } else {
            ss_iter->first -= frontOffset;
        }
#ifdef DEBUG
        if(ss_iter->second > static_cast<unsigned int>(RH_Seq.length())) {
			std::stringstream ss;
			ss<<"Something wrong with front offset! "
            <<"ss iter: "<<ss_iter->first << " \n "
            <<"front offset: " << frontOffset;
			throw crispr::exception(__FILE__,
			                        __LINE__,
			                        __PRETTY_FUNCTION__,
			                        (ss.str()).c_str());
		}
#endif
        ss_iter->second = static_cast<unsigned int>(ss_iter->first + usable_length);
        if (ss_iter->second >= RH_Seq.length()) {
            ss_iter->second = static_cast<unsigned int>(RH_Seq.length() - 1);
        }
    }

    
    // now we check to see if we can find one more DRs on the front or back of this mofo
    // front first
    ss_iter = RH_StartStops.begin();
    if((*ss_iter).first > opts->lowSpacerSize)
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
        
		stringPair sp = smithWaterman(RH_Seq, *DR, &part_s, &part_e, 0, (static_cast<int>((*ss_iter).first) - opts->lowSpacerSize), CRASS_DEF_PARTIAL_SIM_CUT_OFF);
		if(0 != part_e)
		{
			if (part_e - part_s >= CRASS_DEF_MIN_PARTIAL_LENGTH)
			{
				if(((DR->rfind(sp.second) + (sp.second).length()) == DR->length()) && (0 == part_s))
				{
#ifdef DEBUG
					logInfo("adding direct repeat to start",10);
					logInfo(sp.first << " : " << sp.second << " : " << part_s << " : " << part_e,10);
					if(part_e < 0) {
						std::stringstream ss;
						ss<<"Adding negative to SS list! " << part_e;
						throw crispr::exception(__FILE__,
						                        __LINE__,
						                        __PRETTY_FUNCTION__,
						                        (ss.str()).c_str());
					}
					if(part_e > (int)RH_Seq.length()) {
						std::stringstream ss;
						ss <<"SS longer than read: " << part_e;
						throw crispr::exception(__FILE__,
						                        __LINE__,
						                        __PRETTY_FUNCTION__,
						                        (ss.str()).c_str());
					}
#endif
					RH_StartStops.insert(RH_StartStops.begin(), std::make_pair(0, static_cast<unsigned int>(part_e)));
				}
			}
		}
    }
    // then the back
    unsigned int end_dist = static_cast<unsigned int>(RH_Seq.length()) - RH_StartStops.back().second;
    if(end_dist > (unsigned int)(opts->lowSpacerSize))
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
        
		stringPair sp = smithWaterman(RH_Seq,
		                              *DR,
		                              &part_s,
		                              &part_e,
		                              (RH_StartStops.back().second + opts->lowSpacerSize),
		                              (end_dist - opts->lowSpacerSize),
		                              CRASS_DEF_PARTIAL_SIM_CUT_OFF);
		if(0 != part_e)
		{
			if (part_e - part_s >= CRASS_DEF_MIN_PARTIAL_LENGTH)
			{
				if((((int)(RH_Seq.length()) - 1 ) == part_e) && (0 == DR->find(sp.second)))
				{
#ifdef DEBUG
					logInfo("adding partial direct repeat to end",10);
					logInfo(sp.first << " : " << sp.second << " : " << part_s << " : " << part_e,10);
					logInfo((int)sp.first.length() - (int)sp.second.length(),10);
#endif
					// in most cases the right index is returned however
					// if the length of the smith waterman alignment differ the index needs to be corrected
					startStopsAdd(part_s + abs((int)sp.first.length() - (int)sp.second.length()), part_e);
				}
			}
		}
    }
}

std::string crass::ReadHolder::DRLowLexi() {
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    int num_repeats = numRepeats();
    // make sure that tere is 4 elements in the array, if not you can only cut one
    if (num_repeats == 1)
    {
        tmp_dr = repeatStringAt(0);
        rev_comp = reverseComplement(tmp_dr);
    }
    else if (2 == num_repeats)
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // take the second
        if (RH_StartStops.front().first == 0)
        {
            tmp_dr = repeatStringAt(1);
            rev_comp = reverseComplement(tmp_dr);
        }
        
        // take the first
        else if (RH_StartStops.back().second == static_cast<unsigned int>(RH_Seq.length()))
        {
            tmp_dr = repeatStringAt(0);
            rev_comp = reverseComplement(tmp_dr);
        }
        // if they both are then just take whichever is longer
        else
        {
            int lenA = RH_StartStops.front().second - RH_StartStops.front().first;
            int lenB = RH_StartStops.back().second - RH_StartStops.back().first;
            
            if (lenA > lenB)
            {
                tmp_dr = repeatStringAt(0);
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = repeatStringAt(1);
                rev_comp = reverseComplement(tmp_dr);
            }
        }
    }
    // long read more than two repeats
    else
    {
        // take the second
        tmp_dr = repeatStringAt(1);
        rev_comp = reverseComplement(tmp_dr);
        
    }
    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        RH_WasLowLexi = true;
#ifdef DEBUG
        logInfo("DR in low lexi"<<endl<<RH_Seq, 9);
#endif
        return tmp_dr;
    }
    else
    {
        reverseComplementSeq();
        RH_WasLowLexi = false;
#ifdef DEBUG
        logInfo("DR not in low lexi"<<endl<<RH_Seq, 9);
#endif
        return rev_comp;
    }
}

void crass::ReadHolder::reverseComplementSeq() {
    //-----
    // Reverse complement the read and fix the start stops
    //
    
    RH_Seq = reverseComplement(RH_Seq);
	if(RH_Seq.empty()) {
		throw crispr::runtime_exception(__FILE__,
		                                __LINE__,
		                                __PRETTY_FUNCTION__,
		                                "Sequence corrupted during reverse complement!"
		                                );
	}
    reverseStartStops();
    RH_WasLowLexi = !RH_WasLowLexi;
}

void crass::ReadHolder::encode() {
    //-----
    // Encode the sequence using RLE
    // Yo FOOL! ONly call this mofo B4 you make any start stops!
    //
    if(RH_StartStops.size() != 0)
        throw crispr::exception(__FILE__,
                                __LINE__,
                                __PRETTY_FUNCTION__,
                                "Trying to squeeze on non-empty start stops!");
    if (this->RH_isSqueezed)
    {
        return;
    }
    else
    {
        std::stringstream rle, seq;
        rle<<this->RH_Seq[0];
        seq<<this->RH_Seq[0];
		//std::cout<<"seq length: "<<RH_Seq.length()<<std::endl;
		//std::cout<<"orig seq: "<<RH_Seq<<std::endl;
		int length = static_cast<int>(this->RH_Seq.length());
        for (int  i = 1; i < length; i++)
        {
            if (this->RH_Seq[i] == this->RH_Seq[i - 1])
            {
                int count = 0;
                do {
                    count++;
                    i++;
                } while ((i < length) && (this->RH_Seq[i] == this->RH_Seq[i - 1]));
                
				if(i < length) {
                	rle << count << this->RH_Seq[i];
                	seq<<this->RH_Seq[i];
					//std::cout<<"b-index: "<<i<<" base: "<<this->RH_Seq[i]<<std::endl;
				} else {
					rle <<count;
					/*throw crispr::runtime_exception(__FILE__,
                     __LINE__,
                     __PRETTY_FUNCTION__,
                     "Assigning index past end of string");
                     */
				}
            }
            else
            {
                rle << this->RH_Seq[i];
                seq << this->RH_Seq[i];
				//std::cout<<"e-index: "<<i<<" base: "<<this->RH_Seq[i]<<std::endl;
            }
        }
		//std::cout<<"seq: \""<<seq.str().c_str()<<'"'<<std::endl;
		//std::cout<<"rle: \""<<rle.str().c_str()<<'"'<<std::endl;
        this->RH_Seq = seq.str();
        this->RH_Rle = rle.str();
        this->RH_isSqueezed = true;
    }
}

void crass::ReadHolder::decode() {
    //-----
    // Go from RLE to normal
    // Call it anytime. Fixes start stops
    //
    std::string tmp = this->expand(true);
    this->RH_isSqueezed = false;
    this->RH_Seq = tmp;
}


std::string crass::ReadHolder::expand(bool fixStopStarts)
{
    //----
    // Expand the string from RLE and fix stope starts as needed
    //
    
    if (!this->RH_isSqueezed)
    {
        return this->RH_Seq;
    }
    else
    {
        std::stringstream tmp;
        
        int main_index = 0;
        unsigned int new_index = 0;
        unsigned int old_index = 0;
        int stop_index = static_cast<int>(RH_Rle.length());
        int next_ss_index = -1;
        crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
        
        // no point in fixing starts and stops if there are none
        if(RH_StartStops.size() == 0) { fixStopStarts = false; }
        
        if(fixStopStarts)
        {
            next_ss_index = (*ss_iter).first;
        }
        
        while (main_index < stop_index)
        {
            if (isdigit(RH_Rle[main_index]))
            {
                int count = RH_Rle[main_index] - '0';
                new_index += count;
                while (count != 0)
                {
                    tmp << RH_Rle[main_index -1];
                    count--;
                }
            }
            else
            {
                if(next_ss_index == static_cast<int>(old_index))
                {
                    if (ss_iter->first == old_index ) {
                        ss_iter ->first = new_index;
                        next_ss_index = ss_iter->second;
                    } else if (ss_iter->second == old_index) {
                        ss_iter->second = new_index;
                        ss_iter++;
                        if (ss_iter < RH_StartStops.end()) {
                            next_ss_index = ss_iter->first;
                        } else {
                            next_ss_index = -1;
                        }
                    } else {
                        throw crispr::exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, "neither the first or second match the old index");
                    }
                }
                tmp << RH_Rle[main_index];
                old_index++;
                new_index++;
            }
            main_index++;
        }
        return tmp.str();
    }
}

bool crass::ReadHolder::getFirstDR(std::string * retStr)
{
    //-----
    // cut the first DR or return false if it all stuffs up
    //
    RH_LastDREnd = 0;
    return getNextDR(retStr);
}

bool crass::ReadHolder::getNextDR(std::string * retStr)
{
    //-----
    // cut the next DR or return false if it all stuffs up
    //
    // make the iterator point to the start of the next DR
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin() + RH_LastDREnd;
    
    // find out where to start and stop the cuts
    int start_cut = -1;
    int end_cut = -1;
    
    if(ss_iter < RH_StartStops.end())
    {
        start_cut = ss_iter->first;
        end_cut = ss_iter ->second;
    }
    else {
        return false;
    }
    
    // check to see if we made any good of start and end cut
    int dist = end_cut - start_cut;
    if(0 != dist)
    {
        *retStr = RH_Seq.substr(start_cut, dist + 1);
        RH_LastDREnd++;
        return true;
    }
    else
    {
        return false;
    }
}

bool crass::ReadHolder::getFirstSpacer(std::string * retStr)
{
    //-----
    // cut the first Spacer or return false if it all stuffs up
    //
    RH_NextSpacerStart = 0;
    try {
        return getNextSpacer(retStr);
    } catch (crispr::substring_exception& e) {
        std::cerr<<e.what()<<std::endl;
        throw crispr::exception(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Failed to get first spacer");
    }
}

bool crass::ReadHolder::getNextSpacer(std::string * retStr)
{
    //-----
    // cut the next Spacer or return false if it all stuffs up
    //
	
	// first check to see if our index offset makes any sense at all
	if(RH_NextSpacerStart > static_cast<int>(RH_StartStops.numberOfSpacers()))
	{
        //std::stringstream ss;
        //ss << "Next spacer start is greater than length "<<RH_NextSpacerStart<<" > "<<RH_StartStops.size() - 1;
        //throw crispr::exception(__FILE__,
        //                        __LINE__,
        //                        __PRETTY_FUNCTION__,
        //                        ss.str().c_str());
		return false;
	}
	
	// get an iterator into the ss list
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
    
    if(0 == RH_NextSpacerStart)
    {
    	// first run
    	if(0 != ss_iter->first)
    	{
    		// read starts with a spacer
    		// the next spacer starts after the first DR
            try {
                *retStr = RH_Seq.substr(0, ss_iter->first);
            } catch (std::out_of_range& e) {
                throw crispr::substring_exception(e.what(), RH_Seq.c_str(), 0, ss_iter->first, __FILE__, __LINE__, __PRETTY_FUNCTION__);
            }
    		RH_NextSpacerStart = 1;
    	}
    	else
    	{
    		// read starts with a DR
            int start_cut = (*ss_iter).second + 1;
    		ss_iter++;
    		if(ss_iter < RH_StartStops.end())
    		{
                try {
                    *retStr = RH_Seq.substr(start_cut, ss_iter->first - start_cut);
                } catch (std::out_of_range& e) {
                    throw crispr::substring_exception(e.what(), RH_Seq.c_str(), start_cut, (ss_iter->first - start_cut), __FILE__, __LINE__, __PRETTY_FUNCTION__);
                }
            }
    		else
            {
                // only one DR in thie whole guy!
                try {
                    
                    *retStr = RH_Seq.substr(start_cut, RH_Seq.length() - start_cut);
                } catch (std::exception& e) {
                    throw crispr::substring_exception(e.what(), RH_Seq.c_str(), start_cut, (int)(RH_Seq.length() - start_cut), __FILE__, __LINE__, __PRETTY_FUNCTION__);
                    
                }
    		}
            RH_NextSpacerStart = 1;
    	}
        return true;
    }
    else
    {
    	// bump the iterator up!
    	
        ss_iter += RH_NextSpacerStart;
    	// we've been here before
    	if(RH_NextSpacerStart == ((int)(RH_StartStops.size()) - 1))
    	{
    		// last one
            if(ss_iter->second < (RH_Seq.length() - 1))
            {
            	// read ends with a spacer
                try {
                    *retStr = RH_Seq.substr(ss_iter->second + 1);
                    RH_NextSpacerStart++;
                    return true;
                } catch (std::exception& e) {
                    throw crispr::substring_exception(e.what(),
                                                      RH_Seq.c_str(),
                                                      0,
                                                      ss_iter->second,
                                                      __FILE__,
                                                      __LINE__,
                                                      __PRETTY_FUNCTION__);
                }
            }
            else
            {
            	// read ends with a DR. No more spacers to get
#ifdef DEBUG
            	if(ss_iter->second > (RH_Seq.length() - 1))
            	{
                    this->printContents(std::cerr);
                    logError( "ss list out of range; "<<ss_iter->second<< " > "<<RH_Seq.length() - 1);
            	}
#endif
                return false;
            }
    	}
    	else
    	{
    		// middle one
            int start_cut = (*ss_iter).second + 1;
            ss_iter++;
    		int length = ss_iter->first - start_cut;
            try {
                *retStr = RH_Seq.substr(start_cut, length);
    		    RH_NextSpacerStart++;
                return true;
            } catch (std::exception& e) {
                throw crispr::substring_exception(e.what(),
                                                  RH_Seq.c_str(),
                                                  0,
                                                  ss_iter->first,
                                                  __FILE__,
                                                  __LINE__,
                                                  __PRETTY_FUNCTION__);
            }
        }
    }
    // should we get here?
    return false;
}

std::string crass::ReadHolder::splitApart(void)
{
    //-----
    // produce a string of the read split into DR and spacers
    //
    stringstream ss;
    std::string working_str;
    std::string sep_str = " ";
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
    try {
        if(0 == ss_iter->first)
        {
            // start with a DR
            if(getFirstDR(&working_str))
            {
                ss << "DR: " << working_str;
            }
            else
            {
                return "---";
            }
            ss << sep_str;
            
            if(getFirstSpacer(&working_str))
            {
                ss << "SP: " << working_str;
            }
            else
            {
                return "---";
            }
            ss << sep_str;
            
            while(1)
            {
                if(getNextDR(&working_str))
                {
                    ss << "DR: " << working_str;
                }
                else
                {
                    break;
                }
                ss << sep_str;
                if(getNextSpacer(&working_str))
                {
                    ss  << "SP: "<< working_str;
                }
                else
                {
                    break;
                }
                ss << sep_str;
            }
        }
        else
        {
            // start with a spacer
            if(getFirstSpacer(&working_str))
            {
                ss << "SP: " << working_str;
            }
            else
            {
                return "---";
            }
            ss << sep_str;
            
            if(getFirstDR(&working_str))
            {
                ss << "DR: " << working_str;
            }
            else
            {
                return "---";
            }
            ss << sep_str;
            
            // oooohh naughty
            while(1)
            {
                if(getNextSpacer(&working_str))
                {
                    ss << "SP: " << working_str;
                }
                else
                {
                    break;
                }
                ss << sep_str;
                if(getNextDR(&working_str))
                {
                    ss << "DR: " << working_str;
                }
                else
                {
                    break;
                }
                ss << sep_str;
            }
        }
    } catch (crispr::substring_exception& e) {
        std::cerr<<e.what()<<std::endl;
        exit(99);
    }
    
    
    return ss.str();
}

std::string crass::ReadHolder::splitApartSimple(void)
{
    //-----
    // produce a string of the read split into DR and spacers
    // without using the nextDR, nextSpacer functions.
    //
    stringstream ss;
    std::string sep_str = " ";
    unsigned int prev_end = 0;
    
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        if(0 == ss_iter->first)
        {
            // starts with a DR
            ss << "DR: " << RH_Seq.substr(0, ss_iter->second + 1) << sep_str;
            prev_end = ss_iter->second;
        }
        else
        {
            // starts with a spacer
            int length = -1;
            if(0 == prev_end)
            {
                length = ss_iter->first;
            }
            else
            {
                length = ss_iter->first - prev_end  - 1;
                prev_end++;
            }
            ss << "SP: " << RH_Seq.substr(prev_end, length) << sep_str;
            int start = ss_iter->first;
            ss << "DR: " << RH_Seq.substr(start, ss_iter->second - start + 1) << sep_str;
            prev_end = ss_iter->second;
            
        }
        
        if(RH_StartStops.end() == (ss_iter + 1))
        {
            // this is the last one.
            if((RH_Seq.length() - 1) != ss_iter->second)
            {
                // ends on spacer
                ss << "SP: " << RH_Seq.substr(prev_end + 1);
            }
        }
        ss_iter++;
    }
    
    return ss.str();
}

void crass::ReadHolder::printContents(void)
{
    this->printContents(std::cout);    
}
void crass::ReadHolder::printContents(std::ostream& out)
{
    //-----
    // la!
    //
    out <<"Header: "<< RH_Header <<"\n"
    << "LowLexi: " << RH_WasLowLexi << " \n ";
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        out << ss_iter->first << ","<<ss_iter->second;
        ss_iter++;
    }
    out << std::endl;
    out << "Sequence:"<< RH_Seq << std::endl;
    out << "Len: " << RH_Seq.length() << std::endl;
    out << "---------------------------------------------" << std::endl;
    out << "---------------------------------------------" << std::endl;
    
}

void crass::ReadHolder::logContents(int logLevel)
{
    //-----
    // LA!
    //
    stringstream ss;
    ss << RH_Header << " -- " << RH_WasLowLexi << " -- ";
    crispr::RepeatArray<unsigned int>::iterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        ss << ss_iter->first << ","<<ss_iter->second;
        ss_iter++;
    }
    std::string bob;
    ss >> bob;
    logInfo(bob, logLevel);
    logInfo(RH_Seq, logLevel);
}


std::ostream& crass::ReadHolder::print (std::ostream& s)
{
#ifdef OUTPUT_READS_FASTQ
    if (RH_IsFasta)
    {
#endif
        s<<'>'<<RH_Header;
        if (RH_Comment.length() > 0)
        {
            s<<'_'<<RH_Comment;
        }
        s<<std::endl<<RH_Seq;
#ifdef OUTPUT_READS_FASTQ
    }
    else
    {
        s<<'@'<<RH_Header<<std::endl<<RH_Seq<<std::endl<<'+';
        if (RH_Comment.length() > 0)
        {
            s<<RH_Comment<<std::endl;
        }
        s<<RH_Qual;
    }
#endif
    return s;
}

// overloaded operators
std::ostream& operator<< (std::ostream& s,  crass::ReadHolder& c)
{
    return c.print(s);
    
}
