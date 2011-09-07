// File: CrisprNode.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Individual nodes for each kmer!
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
#include <vector>
#include <string>
#include <sstream>
// local includes
#include "CrisprNode.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "GraphDrawingDefines.h"

// boolean value to tell us whether the current node is the first or second in the partner pair
// that is to say when the reads are in lowlexi, which node was seen first in the read
bool CrisprNode::addPartner(CrisprNode * partnerNode, bool nf)
{
    partnerList::iterator part_iter = CN_partnerIDs.begin();
    while (part_iter != CN_partnerIDs.end())
    {
        if (part_iter->CP_partnerNode == partnerNode)
        {
            return false;
        }
        part_iter++;
    }
    
    crisprPartner tmp_cp = { partnerNode, nf };
    CN_partnerIDs.push_back(tmp_cp);
    return true;
}

// print the edges so that the first member of the pair is first
void CrisprNode::printEdges(void)
{
    partnerList::iterator part_iter = CN_partnerIDs.begin();
    while (part_iter != CN_partnerIDs.end())
    {
        // check to see who comes first and print accorndingly
        if (part_iter->CP_nodeFirst) 
        {
            gvEdge(std::cout,part_iter->CP_partnerNode->CN_id, this->CN_id)
        } 
        else 
        {
            gvEdge(std::cout,this->CN_id, part_iter->CP_partnerNode->CN_id)
        }        
        
        part_iter++;
    }
}
