// File: ReadHolder.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Holder of reads. Identified by the various serach algorithms
// Storage class, so stupidly public!
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

#ifndef ReadHolder_h
#define ReadHolder_h

// system includes
#include <iostream>
#include <vector>

// local includes
#include "crass_defines.h"

// typedefs
typedef std::vector<unsigned int> StartStopList;
typedef std::vector<unsigned int>::iterator StartStopListIterator;

class ReadHolder 
{
    public:
        ReadHolder() {}  
        ~ReadHolder() {}  
    
        // members
        std::string RH_Header;                // header for the sequence
        std::string RH_Seq;                   // The DR_lowlexi sequence of this read
        bool RH_WasLowLexi;                   // was the sequence DR_low lexi in the file?
        StartStopList RH_StartStops;          // start stops for DRs, (must be even in length!)
};

#endif //ReadHolder_h