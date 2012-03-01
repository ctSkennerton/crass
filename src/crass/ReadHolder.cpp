// File: ReadHolder.cpp
// Original Author: Michael Imelfort 2011
// Hacked and Extended: Connor Skennerton 2011
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
#include "Exception.h"



// the input must be an even number which will be the start of the repeat
unsigned int ReadHolder::getRepeatAt(unsigned int i)
{
#ifdef DEBUG
    if (i % 2 != 0) 
    {
        logError("Attempting to get a repeat start index from a odd numbered index in RH_StartStops: "<<i);
    }
    if (i > RH_StartStops.size()) 
    {
        logError("Index is greater than the length of the Vector: "<<i);
    }
#endif
    return RH_StartStops[i];
}
std::string ReadHolder::repeatStringAt(unsigned int i)
{
#ifdef DEBUG
	if (i % 2 != 0) 
    {
        logError("Attempting to cut a repeat sequence from a odd numbered index in RH_StartStops: "<<i)
    }
    if (i > RH_StartStops.size()) 
    {
        logError("Index is greater than the length of the Vector: "<<i);
    }
#endif
    return RH_Seq.substr(RH_StartStops[i], RH_StartStops[i + 1] - RH_StartStops[i] + 1);
}

std::string ReadHolder::spacerStringAt(unsigned int i)
{
#ifdef DEBUG
    if (i % 2 != 0) 
    {
        logError("Attempting to cut a spacer sequence from a odd numbered index in RH_StartStops: "<<i)
    }
    if (i > RH_StartStops.size()) 
    {
        logError("Index is greater than the length of the Vector: "<<i);
    }
#endif
    unsigned int curr_spacer_start_index = RH_StartStops[i + 1] + 1;
    unsigned int curr_spacer_end_index = RH_StartStops[i + 2] - 1;
    return RH_Seq.substr(curr_spacer_start_index, (curr_spacer_end_index - curr_spacer_start_index));
}

unsigned int ReadHolder::getAverageSpacerLength()
{
	//-----
	// Used to determine if the string seaching functions are doing a good job
	// eats poo and dies when the read starts / ends with a spacer. So we need to watch that
	//
    unsigned int sum = 0;
    int stored_len = 0;
    unsigned int num_spacers = 0;

    // get the first spacer, but only include it if a DR lies before it
    std::string tmp_string;
    if(getFirstSpacer(&tmp_string))
    {		
        // check to make sure that the read doesn't start on a spacer
        if(*(RH_StartStops.begin()) == 0)
        {
            // starts on a DR
            stored_len = (int)tmp_string.length();
        }
        // else stored_len == 0
        
        // get all the middle spacers
        while(getNextSpacer(&tmp_string))
        {
            num_spacers++;
            sum += stored_len;
            stored_len = (int)tmp_string.length();
        }
        // check to make sure that the read doesn't end on a spacer
        if(RH_StartStops.back() == (RH_Seq.length() - 1))
        {
            // ends on a DR
            num_spacers++;
            sum += stored_len;
        }
        
        if(0 != num_spacers)
            return sum/numSpacers();
    }
    else
    {
        logError("No Spacers!");
    }
    return 0;


    		

}

std::vector<std::string> ReadHolder::getAllSpacerStrings(void)
{
    std::vector<std::string> spacer_vec;
    unsigned int final_index = (getStartStopListSize()/2) + 1;
    for (unsigned int i = 0; i <= final_index; i+=2)
    {
        spacer_vec.push_back(spacerStringAt(i));
    }
    return spacer_vec;
}

std::vector<std::string> ReadHolder::getAllRepeatStrings(void)
{
    std::vector<std::string> repeat_vec;
    unsigned int final_index = getStartStopListSize() - 2;
    for (unsigned int i = 0; i < final_index; i+=2)
    {
        
        repeat_vec.push_back(repeatStringAt(i));
    }
    return repeat_vec;
}

int ReadHolder::averageRepeatLength()
{
    int sum = 0;
    unsigned int final_index = getStartStopListSize() - 2;
    for (unsigned int i = 0; i < final_index; i+=2)
    {
        sum += (int)repeatStringAt(i).length();
    }
    return sum/numRepeats();
}

void ReadHolder::startStopsAdd(unsigned int i, unsigned int j)
{
#ifdef DEBUG	
	if(((int)i < 0) || ((int)j < 0)) { logError("Adding negative to SS list! " << i << " : " << j);	}
	if(i > j) { logError("SS list corrupted! " << i << " : " << j);	}
	if((i > RH_Seq.length()) || (j > RH_Seq.length())) { logError("Too long! " << i << " : " << j);	}
#endif
    this->RH_StartStops.push_back(i);
    if(j >= (unsigned int)getSeqLength())
    {
    	j = (unsigned int)getSeqLength() - 1;
    }
    this->RH_StartStops.push_back(j);
}

void ReadHolder::dropPartials(void)
{
    //-----
    // Drop any partial DRs
    //
    StartStopListIterator r_iter = RH_StartStops.begin();
    if(*r_iter == 0)
    {
        logInfo("\tDropping front partial repeat", 8);
        // this is a partial
        RH_StartStops.erase(r_iter, r_iter+2);
    }
    
    r_iter = RH_StartStops.end() - 1;
    if(*r_iter > (unsigned int)RH_Seq.length() - RH_RepeatLength)
    {
        logInfo("\tDropping end partial repeat", 8);
        // this is a partial
        // TODO This -2 may only need to be -1
        RH_StartStops.erase(r_iter-1, RH_StartStops.end());
    }
}

void ReadHolder::reverseStartStops(void)
{
    //-----
    // the starts and stops were made before the read
    // was potentially rev_comped.
    // So they may need to be fixed
    //
    // we need to fix this mo-fo
    
    StartStopList tmp_ss;
    
    int seq_len = (int)RH_Seq.length();
    int true_start_offset = seq_len - RH_StartStops.back() - 1;
#ifdef DEBUG
    if (true_start_offset < 0) 
    {
        logError("The first direct repeat position is a negative number");
    }
#endif
    StartStopListRIterator ss_iter = RH_StartStops.rbegin();
    
    unsigned int prev_pos_fixed = true_start_offset;
    unsigned int prev_pos_orig = *ss_iter;
   
    while(ss_iter != RH_StartStops.rend())
    {
        unsigned int gap = prev_pos_orig - *ss_iter;
        prev_pos_fixed += gap;
#ifdef DEBUG
        if((int)prev_pos_fixed < 0) { logError("About to add negative value! " << true_start_offset << " : " << prev_pos_orig << " : " << *ss_iter << " : " << gap << " : " << prev_pos_fixed); printContents(); }
#endif
        tmp_ss.push_back(prev_pos_fixed);
        prev_pos_orig = *ss_iter;
        ss_iter++;
    }
    RH_StartStops.clear();
    RH_StartStops.insert(RH_StartStops.begin(), tmp_ss.begin(), tmp_ss.end());
}

void ReadHolder::updateStartStops(int frontOffset, std::string * DR, const options * opts)
{
    //-----
    // Update the start and stops to capture the largest part
    // of the DR
    //
    // Take this opportunity to look for partials at either end of the read
    //
    
    int DR_length = (int)DR->length();
    
    StartStopListIterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {

        int usable_length = DR_length - 1;
        
        // the first guy is the start of the DR
        if(frontOffset >= (int)(*ss_iter))
        {
            // this will be less than 0
            int amount_below_zero = frontOffset - (int)(*ss_iter);
            usable_length = DR_length - amount_below_zero - 1;
            *ss_iter = 0;
        }
        else
        {
            *ss_iter -= frontOffset;
        }
#ifdef DEBUG
        if(*ss_iter > (unsigned int)RH_Seq.length()) { logError("Something wrong with front offset! " << *ss_iter << " : " << frontOffset); }
#endif
        // the second guy is the end of the DR
        ss_iter++;
        *ss_iter = *(ss_iter - 1) + usable_length;
        
        // correct if we have gone beyond the end of the read
        if(*ss_iter >= RH_Seq.length())
        {
            *ss_iter = (unsigned int)RH_Seq.length() - 1;
        }
        ss_iter++;
        
    }

    // now we check to see if we can find one more DRs on the front or back of this mofo
    // front first
    ss_iter = RH_StartStops.begin();
    if((*ss_iter) > opts->lowSpacerSize)
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
        try {
            stringPair sp = smithWaterman(RH_Seq, *DR, &part_s, &part_e, 0, ((int)(*ss_iter) - opts->lowSpacerSize), CRASS_DEF_PARTIAL_SIM_CUT_OFF);
            if(0 != part_e)
            {
                if (part_e - part_s >= CRASS_DEF_MIN_PARTIAL_LENGTH) 
                {
                    if(((DR->rfind(sp.second) + (sp.second).length()) == DR->length()) && (0 == part_s))
                    {
#ifdef DEBUG
                        logInfo("adding direct repeat to start",10);
                        logInfo(sp.first << " : " << sp.second << " : " << part_s << " : " << part_e,10);
                    	if(part_e < 0) { logError("Adding negative to SS list! " << part_e); }
                    	if(part_e > (int)RH_Seq.length()) { logError("SS longer than read: " << part_e); }
#endif
                        std::reverse(RH_StartStops.begin(), RH_StartStops.end());
                        RH_StartStops.push_back(part_e);
                        RH_StartStops.push_back(0);
                        std::reverse(RH_StartStops.begin(), RH_StartStops.end());
                    }
                }
            }
        } catch (crispr::exception& e) {
            std::cerr<<e.what()<<std::endl;
            exit(148);
        }
    }
    // then the back
    unsigned int end_dist = (unsigned int)RH_Seq.length() - RH_StartStops.back();
    if(end_dist > (unsigned int)(opts->lowSpacerSize))
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
        try {
            stringPair sp = smithWaterman(RH_Seq, *DR, &part_s, &part_e, (RH_StartStops.back() + opts->lowSpacerSize), (end_dist - opts->lowSpacerSize), CRASS_DEF_PARTIAL_SIM_CUT_OFF);
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
        } catch (crispr::exception& e) {
            std::cerr<<e.what()<<std::endl;
            exit(148);
        }
    }
}

std::string ReadHolder::DRLowLexi(void)
{
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
        if (RH_StartStops.front() == 0)
        {
            tmp_dr = repeatStringAt(2);
            rev_comp = reverseComplement(tmp_dr);
        }
        
        // take the first
        else if (RH_StartStops.back() == (unsigned int)RH_Seq.length())
        {
            tmp_dr = repeatStringAt(0);
            rev_comp = reverseComplement(tmp_dr);
        }
        // if they both are then just take whichever is longer
        else
        {
            int lenA = RH_StartStops.at(1) - RH_StartStops.at(0);
            int lenB = RH_StartStops.at(3) - RH_StartStops.at(2);
            
            if (lenA > lenB)
            {
                tmp_dr = repeatStringAt(0);
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = repeatStringAt(2);
                rev_comp = reverseComplement(tmp_dr);
            }
        }
    }
    // long read more than two repeats
    else
    {
        // take the second
        tmp_dr = repeatStringAt(2);
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

void ReadHolder::reverseComplementSeq(void)
{
    //-----
    // Reverse complement the read and fix the start stops
    // 
    RH_Seq = reverseComplement(RH_Seq);
    reverseStartStops();
    RH_WasLowLexi = !RH_WasLowLexi;
}

// simple run length encoding
void ReadHolder::encode(void)
{
    //-----
    // Encode the sequence using RLE
    // Yo FOOL! ONly call this mofo B4 you make any start stops!
    //
    if(RH_StartStops.size() != 0)
        logError("Trying to squeeze on non-empty start stops!");
    if (this->RH_isSqueezed) 
    {
        return;
    } 
    else 
    {
        std::stringstream rle, seq;
        rle<<this->RH_Seq[0];
        seq<<this->RH_Seq[0];
        for (int  i = 1; i < (int)(this->RH_Seq.length()); i++) 
        {
            if (this->RH_Seq[i] == this->RH_Seq[i - 1]) 
            {
                int count = 0;
                do {
                    count++;
                    i++;
                } while (this->RH_Seq[i] == this->RH_Seq[i - 1]); 
                
                rle << count << this->RH_Seq[i];
                seq<<this->RH_Seq[i];
            }
            else
            {
                rle << this->RH_Seq[i];
                seq << this->RH_Seq[i];
            }
        }
        this->RH_Seq = seq.str();
        this->RH_Rle = rle.str();
        this->RH_isSqueezed = true;
    }
}

void ReadHolder::decode(void)
{
    //-----
    // Go from RLE to normal
    // Call it anytime. Fixes start stops
    //
    std::string tmp = this->expand(true);
    this->RH_isSqueezed = false;
    this->RH_Seq = tmp;
}


std::string ReadHolder::expand(void)
{
    //----
    // Expand the string from RLE but don't fix stop starts.
    //
    return expand(false);
}

std::string ReadHolder::expand(bool fixStopStarts)
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
        int new_index = 0;
        int old_index = 0;
        int stop_index = (int)RH_Rle.length();
        int next_ss_index = -1;
        StartStopListIterator ss_iter = RH_StartStops.begin();
        
        // no point in fixing starts and stops if there are none
        if(RH_StartStops.size() == 0) { fixStopStarts = false; }
        
        if(fixStopStarts)
        {
            next_ss_index = *ss_iter;
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
                if(next_ss_index == old_index)
                {
                    *ss_iter = new_index;
                    ss_iter++;
                    if(ss_iter != RH_StartStops.end())
                    {
                        next_ss_index = *ss_iter;
                    }
                    else
                    {
                        next_ss_index = -1;
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

// cut DRs and Specers

bool ReadHolder::getFirstDR(std::string * retStr)
{
    //-----
    // cut the first DR or return false if it all stuffs up
    //
    RH_LastDREnd = 0;
    return getNextDR(retStr);
}

bool ReadHolder::getNextDR(std::string * retStr)
{
    //-----
    // cut the next DR or return false if it all stuffs up
    //
    // make the iterator point to the start of the next DR
    StartStopListIterator ss_iter = RH_StartStops.begin() + RH_LastDREnd;
    
    // find out where to start and stop the cuts
    int start_cut = -1;
    int end_cut = -1;
    
    if(ss_iter < RH_StartStops.end())
    {
        start_cut = *ss_iter;
    }
    else
        return false;
    ss_iter++;
    if(ss_iter < RH_StartStops.end())
    {
        end_cut = *ss_iter;
    }
    else
    {
        return false;
    }
    
    // check to see if we made any good of start and end cut
    int dist = end_cut - start_cut;
    if(0 != dist)
    {
        *retStr = RH_Seq.substr(start_cut, dist + 1);
        RH_LastDREnd+=2;
        return true;
    }
    else
    {
        return false;
    }
}

bool ReadHolder::getFirstSpacer(std::string * retStr)
{
    //-----
    // cut the first Spacer or return false if it all stuffs up
    //
    RH_NextSpacerStart = 0;
    try {
        return getNextSpacer(retStr);
    } catch (crispr::substring_exception& e) {
        std::cerr<<e.what()<<std::endl;
        exit(99);
    }
}

bool ReadHolder::getNextSpacer(std::string * retStr)
{
    //-----
    // cut the next Spacer or return false if it all stuffs up
    //
	
	// first check to see if our index offset makes any sense at all
	if(RH_NextSpacerStart > ((int)(RH_StartStops.size()) - 1))
	{
		return false;
	}
	
	// get an iterator into the ss list
    StartStopListIterator ss_iter = RH_StartStops.begin();
    
    if(0 == RH_NextSpacerStart)
    {
    	// first run
    	if(0 != *ss_iter)
    	{
    		// read starts with a spacer
    		// the next spacer starts after the first DR
            try {
                *retStr = RH_Seq.substr(0, *ss_iter);
            } catch (std::out_of_range& e) {
                throw crispr::substring_exception(e.what(), RH_Seq.c_str(), 0, *ss_iter, __FILE__, __LINE__, __PRETTY_FUNCTION__);
            }
    		RH_NextSpacerStart = 1;
    	}
    	else
    	{
    		// read starts with a DR
            ss_iter++;
            int start_cut = (*ss_iter) + 1;
    		ss_iter++;
    		if(ss_iter < RH_StartStops.end())
    		{
                try {
                    *retStr = RH_Seq.substr(start_cut, *ss_iter - start_cut);
                } catch (std::out_of_range& e) {
                    throw crispr::substring_exception(e.what(), RH_Seq.c_str(), start_cut, (*ss_iter - start_cut), __FILE__, __LINE__, __PRETTY_FUNCTION__);
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
            RH_NextSpacerStart = 3;
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
            if(*ss_iter < (RH_Seq.length() - 1))
            {
            	// read ends with a spacer
                try {
                    *retStr = RH_Seq.substr(*ss_iter + 1);
                    RH_NextSpacerStart+=2;
                    return true;
                } catch (std::exception& e) {
                    throw crispr::substring_exception(e.what(), RH_Seq.c_str(), 0, *ss_iter, __FILE__, __LINE__, __PRETTY_FUNCTION__);
                }
            }
            else
            {
            	// read ends with a DR. No more spacers to get
#ifdef DEBUG
            	if(*ss_iter > (RH_Seq.length() - 1))
            	{
            		logError("ss list out of range");
            	}
#endif
                return false;    		
            }
    	}
    	else
    	{
    		// middle one
            int start_cut = (*ss_iter) + 1;
            ss_iter++;
    		int length = *ss_iter - start_cut;
            try {
                *retStr = RH_Seq.substr(start_cut, length);
    		    RH_NextSpacerStart += 2;
                return true;
            } catch (std::exception& e) {
                throw crispr::substring_exception(e.what(), RH_Seq.c_str(), 0, *ss_iter, __FILE__, __LINE__, __PRETTY_FUNCTION__);
            }
        }
    }    
    // should we get here?
    return false;
}



std::string ReadHolder::toStringInColumns(void)
{
    //-----
    // produce a string of the read split into DR and spacers
    //

    std::stringstream str;
    
    std::string repeat, spacer, prev_spacer;
    repeat = spacer = prev_spacer = "";
    
    
    str << "POSITION\tREPEAT\t\t\t\tSPACER"<<std::endl;
    
    str <<"--------\t";
    
    for (int y = 0; y <  RH_RepeatLength; y++) str <<"-";
    str<<"\t";
    
    for (unsigned int z = 0; z < getAverageSpacerLength(); z++) str << "-";
    str<<std::endl;
    
    
    //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
    for (unsigned int m = 0; m < getStartStopListSize(); m+=2)
    {   //repeat = getRepeat(m);
        str << RH_StartStops[m] + 1 << "\t\t" << repeatStringAt((unsigned int)m) << "\t";
        
        // print spacer
        // because there are no spacers after the last repeat, we stop early (m < crisprIndexVector.size() - 1)
        if (m < numSpacers())
        {   prev_spacer = spacer;
            spacer = spacerStringAt(m);
            str << spacer;
            
            str <<"\t[ " << repeatStringAt(m).length() << ", " << spacerStringAt(m).length() << " ]";
            str <<std::endl;
            
        }
    }
    
    
    str <<std::endl<<"--------\t";
    
    for (int x = 0; x < RH_RepeatLength; x++)
    {
        str << "-";
    }
    str <<"\t";
    
    for (unsigned int z = 0; z <  getAverageSpacerLength(); z++)
    {
        str << "-";
    }
    str <<std::endl;
    
    
    return str.str();
    
}


std::string ReadHolder::splitApart(void)
{
    //-----
    // produce a string of the read split into DR and spacers
    //
    stringstream ss;
    std::string working_str;
    std::string sep_str = " ";
    StartStopListIterator ss_iter = RH_StartStops.begin();
    try {
        if(0 == *ss_iter)
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

std::string ReadHolder::splitApartSimple(void)
{
    //-----
    // produce a string of the read split into DR and spacers
    // without using the nextDR, nextSpacer functions.
    //
    stringstream ss;
    std::string sep_str = " ";
    unsigned int prev_end = 0;
    
    StartStopListIterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        if(0 == *ss_iter)
        {
            // starts with a DR
            ss_iter++;
            ss << "DR: " << RH_Seq.substr(0, *ss_iter + 1) << sep_str;
            prev_end = *ss_iter;
        }
        else
        {
            // starts with a spacer
            int length = -1;
            if(0 == prev_end)
            {
                length = *ss_iter;
            }
            else
            {
                length = *ss_iter - prev_end  - 1;
                prev_end++;
            }
            ss << "SP: " << RH_Seq.substr(prev_end, length) << sep_str;
            int start = *ss_iter;
            ss_iter++;
            ss << "DR: " << RH_Seq.substr(start, *ss_iter - start + 1) << sep_str;
            prev_end = *ss_iter;
            
        }
        
        if(RH_StartStops.end() == (ss_iter + 1))
        {
            // this is the last one.
            if((RH_Seq.length() - 1) != *ss_iter)
            {
                // ends on spacer
                ss << "SP: " << RH_Seq.substr(prev_end + 1);
            }
        }
        ss_iter++;
    }
    
    return ss.str();
}

void ReadHolder::printContents(void)
{
    //-----
    // la!
    //
    std::cout << RH_Header << " -- " << RH_WasLowLexi << " -- " << std::flush;
    StartStopListIterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        std::cout << *ss_iter << ",";
        ss_iter++;
    }
    std::cout << std::endl;
    std::cout << RH_Seq << std::endl;
    std::cout << "Len: " << RH_Seq.length() << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    
}

void ReadHolder::logContents(int logLevel)
{
    //-----
    // LA!
    //
    stringstream ss;
    ss << RH_Header << " -- " << RH_WasLowLexi << " -- ";
    StartStopListIterator ss_iter = RH_StartStops.begin();
    while(ss_iter != RH_StartStops.end())
    {
        ss << *ss_iter << ",";
        ss_iter++;
    }
    std::string bob;
    ss >> bob;
    logInfo(bob, logLevel);
    logInfo(RH_Seq, logLevel);
}


std::ostream& ReadHolder::print (std::ostream& s)
{
    if (RH_IsFasta) 
    {
        s<<'>'<<RH_Header;
        if (RH_Comment.length() > 0) 
        {
            s<<'_'<<RH_Comment;
        }
        s<<std::endl<<RH_Seq;
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
    return s;
}

// overloaded operators 
std::ostream& operator<< (std::ostream& s,  ReadHolder& c)
{
    return c.print(s);
    
}
