// File: ReadHolder.cpp
// Original Author: Michael Imelfort 2011
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

// local includes
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "SmithWaterman.h"
#include "LoggerSimp.h"


std::string ReadHolder::seq(void)
{
    if (this->RH_isSqueezed) 
    {
        std::stringstream tmp;
        std::string::iterator str_iter = this->RH_Seq.begin();
        
        while (str_iter != this->RH_Seq.end()) 
        {
            if (!isdigit(*str_iter)) 
            {
                tmp<<*str_iter;
            }
            str_iter++;
        }
        return tmp.str();
    }    
    else
    {
        return this->RH_Seq;
    }
}

void ReadHolder::add(int i)
{
    this->RH_StartStops.push_back(i);
}

void ReadHolder::add(int i, int j)
{
    this->RH_StartStops.push_back(i);
    this->RH_StartStops.push_back(j);
}

int ReadHolder::at(int i)
{
    return this->RH_StartStops.at(i);
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
    
    StartStopListRIterator ss_iter = RH_StartStops.rbegin();

    unsigned int prev_pos_fixed = true_start_offset;
    unsigned int prev_pos_orig = *ss_iter;
    
    while(ss_iter != RH_StartStops.rend())
    {
        unsigned int gap = prev_pos_orig - *ss_iter;
        prev_pos_fixed += gap;
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
        if(frontOffset >= *ss_iter)
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
    if((int)(*ss_iter) > opts->lowSpacerSize)
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
//        std::cout << "Check start: " << *ss_iter << " : " << opts->lowSpacerSize << std::endl;
        stringPair sp = smithWaterman(RH_Seq, *DR, &part_s, &part_e, 0, ((int)(*ss_iter) - opts->lowSpacerSize), CRASS_DEF_PARTIAL_SIM_CUT_OFF);
        if(0 != part_e)
        {
            if(((DR->rfind(sp.second) + (sp.second).length()) == DR->length()) && (0 == part_s))
            {
                //std::cout << sp.first << " : " << sp.second << " : " << part_s << " : " << part_e << std::endl;
                std::reverse(RH_StartStops.begin(), RH_StartStops.end());
                RH_StartStops.push_back(part_e);
                RH_StartStops.push_back(0);
                std::reverse(RH_StartStops.begin(), RH_StartStops.end());
            }
        }
    }
    
    // then the back
    unsigned int end_dist = (unsigned int)RH_Seq.length() - RH_StartStops.back();
    if(end_dist > opts->lowSpacerSize)
    {
        // we should look for a DR here
        int part_s, part_e;
        part_s = part_e = 0;
//        std::cout << "Check end: " << end_dist << " : " << opts->lowSpacerSize << std::endl;
        stringPair sp = smithWaterman(RH_Seq, *DR, &part_s, &part_e, (RH_StartStops.back() + opts->lowSpacerSize), (end_dist - opts->lowSpacerSize), CRASS_DEF_PARTIAL_SIM_CUT_OFF);
        if(0 != part_e)
        {
            if(((RH_Seq.length() - 1 ) == part_e) && (0 == DR->find(sp.second)))
            {
                //std::cout << sp.first << " : " << sp.second << " : " << part_s << " : " << part_e << std::endl;
                RH_StartStops.push_back(part_s);
                RH_StartStops.push_back(part_e);
            }
        }
    }
}

// cut DRs and Specers

bool ReadHolder::getFirstDR(std::string * retStr)
{
    //-----
    // cut the first DR or return false if it all stuffs up
    //
    mLastDREnd = 0;
    return getNextDR(retStr);
}

bool ReadHolder::getNextDR(std::string * retStr)
{
    //-----
    // cut the next DR or return false if it all stuffs up
    //
    // make the iterator point to the start of the next DR
    StartStopListIterator ss_iter = RH_StartStops.begin() + mLastDREnd;
    
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
        *retStr = RH_Seq.substr(start_cut, dist+1);
        mLastDREnd+=2;
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
    mLastSpacerEnd = 0;
    return getNextSpacer(retStr);
}

bool ReadHolder::getNextSpacer(std::string * retStr)
{
    //-----
    // cut the next Spacer or return false if it all stuffs up
    //
    StartStopListIterator ss_iter = RH_StartStops.begin();

    // find out where to start and stop the cuts
    int start_cut = -1;
    int end_cut = -1;

    if(mLastSpacerEnd == (RH_StartStops.size() - 1))
    {
        ss_iter += mLastSpacerEnd;
        if(RH_Seq.length() != ((*ss_iter) + 1))
        {
            mLastSpacerEnd++;
            *retStr = RH_Seq.substr(*ss_iter + 1);
            return true;
        }
        return false;
    }

    // if the DR starts at 0, 
    if(0 == mLastSpacerEnd)
    {
        // first run
        if(0 != *ss_iter)
        {
            // cut the front before the first DR
            start_cut = 0;
            end_cut = *ss_iter;
        }
        mLastSpacerEnd = 1;
    }
    
    if(-1 == start_cut)
    {
        // we didn't set it above
        ss_iter = RH_StartStops.begin() + mLastSpacerEnd;
    
        if(ss_iter < RH_StartStops.end())
        {
            start_cut = *ss_iter;
        }
        else
        {
            return false;
        }
        ss_iter++;
        if(ss_iter < RH_StartStops.end())
        {
            end_cut = *ss_iter;
        }
        else
        {
            return false;
        }

        mLastSpacerEnd += 2; 
    }
        
    // check to see if we made any good of start and end cut
    if(0 != start_cut)
    {
        start_cut++;
    }
    if(end_cut == RH_Seq.length() - 1)
    {
        *retStr = RH_Seq.substr(start_cut);
    }
    else
    {
        *retStr = RH_Seq.substr(start_cut, end_cut - start_cut);
    }
    return true;
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

void ReadHolder::encode(void)
{
    std::string tmp = this->squeeze();
    this->RH_isSqueezed = true;
    this->RH_Seq = tmp;
}

void ReadHolder::decode(void)
{
    std::string tmp = this->expand();
    this->RH_isSqueezed = false;
    this->RH_Seq = tmp;
}

// simple run length encoding
std::string ReadHolder::squeeze(void)
{
    if (this->RH_isSqueezed) 
    {
        return this->RH_Seq;
    } 
    else 
    {
        std::stringstream tmp;
        // make sure that the first comparison is to something that is never in DNA
        char previous_base = 'Z';
        
        for (int  i = 0; i < this->RH_Seq.length(); i++) 
        {
            int count = 0;
            if (this->RH_Seq[i] == previous_base) 
            {
                while (this->RH_Seq[i] == previous_base) 
                {
                    count++;
                    i++;
                }
                tmp << count << this->RH_Seq[i];
            }
            else
            {
                tmp << this->RH_Seq[i];
            }

            previous_base = this->RH_Seq[i];
        }
        return tmp.str();
    }
}

std::string ReadHolder::expand(void)
{
    if (!this->RH_isSqueezed) 
    {
        return this->RH_Seq;
    } 
    else 
    {
        std::stringstream tmp;
        std::string::iterator str_iter = this->RH_Seq.begin();
        
        while (str_iter != this->RH_Seq.end()) 
        {
            if (isdigit(*str_iter)) 
            {
                int count = *str_iter - '0';
                while (count != 0) 
                {
                    tmp<<*(str_iter - 1);
                    count--;
                }
            }
            else
            {
                tmp<<*str_iter;
            }
            str_iter++;
        }
        return tmp.str();
    }
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
