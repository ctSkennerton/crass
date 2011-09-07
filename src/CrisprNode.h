// File: CrisprNode.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Header file for the NodeManager
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

#ifndef CrisprNode_h
    #define CrisprNode_h

// system includes
#include <vector>
#include <string>

// local includes
#include "crassDefines.h"
#include "StringCheck.h"

class CrisprNode;

// struct of edges
// bool is true if the partner lies before the current node
typedef struct 
{
    CrisprNode * CP_partnerNode;
    bool CP_nodeFirst;
} crisprPartner;

typedef std::vector<crisprPartner> partnerList;

class CrisprNode 
{
    public:
        //constructor
        CrisprNode(StringToken id)
        {
            CN_id = id;
        }
        //destructor
        ~CrisprNode(){}
        
        // return success if the partner has been added
        // these are basically the edges of the graph
        bool addPartner(CrisprNode * parterNode, bool nf);
        
        partnerList * getPartners(void)
        {
            return &CN_partnerIDs;
        }    
        
    void printEdges(void);
    
    private:
        // id of the kmer of the cripsr node
        StringToken CN_id;
        // the ID of the partner of the current node in the lookup table
        partnerList CN_partnerIDs;

        
};

#endif //CrisprNode_h