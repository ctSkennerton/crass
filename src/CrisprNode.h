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
#include <fstream>

// local includes
#include "crassDefines.h"
#include "StringCheck.h"
#include "Rainbow.h"

class CrisprNode;

// Enum to let us know if the node is a "first" node in a spacer pair
enum EDGE_TYPE {
    CN_EDGE_BACKWARD,
    CN_EDGE_FORWARD,
    CN_EDGE_JUMPING_F,
    CN_EDGE_JUMPING_B
};

// a list of edges ( we use a map to make lookups faster )
// The bool teels us if the edge is active (ie, if the joining node is still attached / in use)
typedef std::map<CrisprNode *, bool> edgeList;
typedef std::map<CrisprNode *, bool>::iterator edgeListIterator;

class CrisprNode 
{
    public:
        //constructor
        CrisprNode(void)
        {
            CN_id = 0;
            CN_Attached = true;
            CN_InnerRank_F = 0;
            CN_InnerRank_B = 0;
            CN_JumpingRank_F = 0;
            CN_JumpingRank_B = 0;
            CN_Coverage = 0;
            CN_IsForward = true;
        }

        CrisprNode(StringToken id)
        {
            CN_id = id;
            CN_Attached = true;                                         // by default, a node is attached unless actually detached by the user
            CN_InnerRank_F = 0;
            CN_InnerRank_B = 0;
            CN_JumpingRank_F = 0;
            CN_JumpingRank_B = 0;
            CN_Coverage = 1;
            CN_IsForward = true;
        }
        
        //destructor
        ~CrisprNode(){}
        
        //
        // Generic get and set
        //
        inline StringToken getID(void) { return CN_id; }
        inline bool isForward(void) { return CN_IsForward; }
        inline void setForward(bool forward) { CN_IsForward = forward; }
        inline int getCoverage() {return CN_Coverage;}
        
        //
        // Edge level functions
        //
        bool addEdge(CrisprNode * parterNode, EDGE_TYPE type);          // return success if the partner has been added
        edgeList * getEdges(EDGE_TYPE type);                            // get edges of a particular type
        
        //
        // Node level functions
        //
        inline void detachNode(void) { setAttach(false); }              // detach this node
        inline void reattachNode(void) { setAttach(true); }             // re-attach this node
        inline bool isAttached(void) { return CN_Attached; }            // der...
        int getRank(EDGE_TYPE type);                                    // return the rank of the node
        inline void incrementCount(void) { CN_Coverage++; }             // Increment the coverage
        
        //
        // File IO / printing
        //

        void printEdges(std::ostream &dataOut, bool showDetached, bool printBackEdges);    

        std::string sayEdgeTypeLikeAHuman(EDGE_TYPE type);
    
    private:
    
        void setAttach(bool attachState);                               // set the attach state of the node
        
        // id of the kmer of the cripsr node
        StringToken CN_id;
        
        //
        // We need different edge lists to store the variety of edges we may encounter, observe...
        //
        //                    NODE 1       NODE 2                 NODE 3       NODE 4                 NODE 5       NODE 6
        //  ... DRDRDRDRDR | SP_start ---- SP_end | DRDRDRDRDR | SP_start ---- SP_end | DRDRDRDRDR | SP_start ---- SP_end | DRDRDRDRDR ...
        //
        // Gives edge types: (B = CN_EDGE_BACKWARD, F = CN_EDGE_FORWARD, JF = JUMPING CN_EDGE_FORWARD, JB = JUMPING CN_EDGE_BACKWARD, X = no egde)
        //
        // ---------------------------------------------------------------
        //         | NODE 1 | NODE 2 | NODE 3 | NODE 4 | NODE 5 | NODE 6 |
        // ---------------------------------------------------------------
        //  NODE 1 |   X    |   F    |   X    |   X    |   X    |   X    |
        //  NODE 2 |   B    |   X    |   JF   |   X    |   X    |   X    |
        //  NODE 3 |   X    |   JB   |   X    |   F    |   X    |   X    |
        //  NODE 4 |   X    |   X    |   B    |   X    |   JF   |   X    |
        //  NODE 5 |   X    |   X    |   X    |   JB   |   X    |   F    |
        //  NODE 6 |   X    |   X    |   X    |   X    |   B    |   X    |
        // ---------------------------------------------------------------
        //
        edgeList CN_ForwardEdges;
        edgeList CN_BackwardEdges;
        edgeList CN_JumpingForwardEdges;
        edgeList CN_JumpingBackwardEdges;

        // We need multiple classes of RANK
        int CN_InnerRank_F;
        int CN_InnerRank_B;
        int CN_JumpingRank_F;
        int CN_JumpingRank_B;

        // We need to mark whether nodes are atached or detached
        bool CN_Attached;
        
        // how many times have we seen this guy?
        int CN_Coverage;
        
        // is this a forward facing node?
        bool CN_IsForward;
};

#endif //CrisprNode_h
