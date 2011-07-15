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

// local includes
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "LoggerSimp.h"

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
    int seq_len = RH_Seq.length();
    int true_start_offset = seq_len - RH_StartStops.back() - 1;
    
    StartStopListIterator ss_iter = RH_StartStops.begin();
    int prev_pos = *ss_iter;
    int pos_holder = true_start_offset;
    int rh_ss_pos = 0;
    ss_iter++;
    while(ss_iter != RH_StartStops.end())
    {
        RH_StartStops[rh_ss_pos] = pos_holder;
        pos_holder = *ss_iter - prev_pos + pos_holder;
        prev_pos = *ss_iter;
        rh_ss_pos++;
        ss_iter++;
    }
    RH_StartStops[rh_ss_pos] = pos_holder;
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
