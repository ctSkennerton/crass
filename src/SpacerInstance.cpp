// File: SpacerInstance.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// Implementation of SpacerInstance functions
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
#include <string>

// local includes
#include "SpacerInstance.h"
#include "CrisprNode.h"
#include "StringCheck.h"

SpacerInstance::SpacerInstance(StringToken spacerID) 
{
    //-----
    // constructor
    //
    SI_SpacerSeqID = spacerID;
    SI_LeadingNode = NULL;
    SI_LastNode = NULL;
    SI_InstanceCount = 0;
    mSpacerRank = 0;
    mContigID = 0;
    mAttached = false;
}

SpacerInstance::SpacerInstance(StringToken spacerID, CrisprNode * leadingNode, CrisprNode * lastNode)
{
    SI_SpacerSeqID = spacerID;
    SI_LeadingNode = leadingNode;
    SI_LastNode = lastNode;
    SI_InstanceCount  = 1;
    mSpacerRank = 0;
    mContigID = 0;
    mAttached = false;
}




