// File: SpacerInstance.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Spacers! Lots of them!
// Each found spacer *could* be different. Each unique instance should be saved
// so that graph cleaning can be doned.
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

#ifndef SpacerInstance_h
    #define SpacerInstance_h

// system includes
#include <iostream>

// local includes
#include "crass_defines.h"
#include "CrisprNode.h"

class SpacerInstance {
    public:
        SpacerInstance (std::string spacerSeq);
        ~SpacerInstance () {}
        
    private:
        std::string mSpacerSeq;                 // the sequence of this spacer
        CrisprNode * mLeadingNode;              // the first node of this spacer
        CrisprNode * mLastNode;                 // the last node
        unsigned int mInstanceCount;            // the number of times this exact instance has been seen
};


#endif //SpacerInstance_h