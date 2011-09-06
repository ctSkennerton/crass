// File: NodeManager.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Implementation of NodeManager functions
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
#include <vector>
#include <string>

// local includes
#include "NodeManager.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "CrisprNode.h"
#include "SpacerInstance.h"
#include "libcrispr.h"
#include "StringCheck.h"
#include "ReadHolder.h"

NodeManager::NodeManager(std::string drSeq)
{
    //-----
    // constructor
    //
    mDirectRepeatSequence = drSeq;
}

NodeManager::~NodeManager(void)
{
    //-----
    // destructor
    //
    
    // clean up al lthe cripsr nodes
    NodeListIterator node_iter = mNodes.begin();
    while(node_iter != mNodes.end())
    {
        if(NULL != *node_iter)
        {
            delete *node_iter;
            *node_iter = NULL;
        }
        node_iter++;
    }
    mNodes.clear();
    
    SpacerListIterator spacer_iter = mSpacers.begin();
    while(spacer_iter != mSpacers.end())
    {
        if(NULL != *spacer_iter)
        {
            delete *spacer_iter;
            *spacer_iter = NULL;
        }
        spacer_iter++;
    }
    mSpacers.clear();
}

bool NodeManager::addReadholder(ReadHolder * RH)
{
    //-----
    // add a readholder to this mofo and fix all the start stops
    //
    mReadList.push_back(RH);
    return true;
}































