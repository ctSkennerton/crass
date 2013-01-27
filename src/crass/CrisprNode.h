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
#include "libcrispr.h"
#include "ReadHolder.h"

class CrisprNode;

// Enum to let us know if the node is a "first" node in a spacer pair
enum EDGE_TYPE {
    CN_EDGE_BACKWARD,
    CN_EDGE_FORWARD,
    CN_EDGE_JUMPING_F,
    CN_EDGE_JUMPING_B,
    CN_EDGE_ERROR
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
            mid = 0;
            mAttached = true;
            mInnerRank_F = 0;
            mInnerRank_B = 0;
            mJumpingRank_F = 0;
            mJumpingRank_B = 0;
            mCoverage = 0;
            mIsForward = true;
        }

        CrisprNode(StringToken id)
        {
            mid = id;
            mAttached = true;                                         // by default, a node is attached unless actually detached by the user
            mInnerRank_F = 0;
            mInnerRank_B = 0;
            mJumpingRank_F = 0;
            mJumpingRank_B = 0;
            mCoverage = 1;
            mIsForward = true;
        }
        
        //destructor
        ~CrisprNode(){}
        
        //
        // Generic get and set
        //
        inline StringToken getID(void) { return mid; }
        inline bool isForward(void) { return mIsForward; }
        inline void setForward(bool forward) { mIsForward = forward; }
        inline int getCoverage() {return mCoverage;}
        int getDiscountedCoverage(void);
        inline void addReadHeader(StringToken readHeader) { mReadHeaders.push_back(readHeader); }
        inline void addReadHolder(crass::ReadHolder * RH) { mReadHolders.push_back(RH); }
        inline std::vector<StringToken> * getReadHeaders(void) { return &mReadHeaders; }
        inline ReadList * getReadHolders(void) { return &mReadHolders; }
        
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
        inline bool isAttached(void) { return mAttached; }            // der...
        void setAsDetached(void) { mAttached = false; }					// DO NOT CALL THIS OUTSIDE OF THE ATTACH FUNCTION!
        int getRank(EDGE_TYPE type);                                    // return the rank of the node
        void updateRank(bool attachState, EDGE_TYPE type);				// increment or decrement the rank of this type
        inline void incrementCount(void) { mCoverage++; }             // Increment the coverage
        int getTotalRank(void) { return getRank(CN_EDGE_BACKWARD) + getRank(CN_EDGE_FORWARD) + getRank(CN_EDGE_JUMPING_F) + getRank(CN_EDGE_JUMPING_B); }
        int getJumpingRank(void) { return getRank(CN_EDGE_JUMPING_F) + getRank(CN_EDGE_JUMPING_B); }
        int getInnerRank(void) { return getRank(CN_EDGE_BACKWARD) + getRank(CN_EDGE_FORWARD); }
        
        //
        // File IO / printing
        //

        void printEdges(std::ostream &dataOut, StringCheck * ST, std::string label, bool showDetached, bool printBackEdges, bool longDesc);    
        std::vector<std::string> getReadHeaders(StringCheck * ST);
        std::string sayEdgeTypeLikeAHuman(EDGE_TYPE type);
    std::vector<StringToken>::iterator beginHeaders(void) {return mReadHeaders.begin();}
    std::vector<StringToken>::iterator endHeaders(void) {return mReadHeaders.end();}

    private:
    
        void setAttach(bool attachState);                               // set the attach state of the node
        void setEdgeAttachState(edgeList * currentList, bool attachState, EDGE_TYPE currentType);
        void calculateReadCoverage(edgeList * currentList, std::map<StringToken, int>& countingMap);
    void printEdgesForList(edgeList * currentList,
                           std::ostream &dataOut, 
                           StringCheck * ST,
                           std::string label, 
                           bool showDetached, 
                           bool longDesc);        
        // id of the kmer of the cripsr node
        StringToken mid;
        
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
        edgeList mForwardEdges;
        edgeList mBackwardEdges;
        edgeList mJumpingForwardEdges;
        edgeList mJumpingBackwardEdges;

        // We need multiple classes of RANK
        int mInnerRank_F;
        int mInnerRank_B;
        int mJumpingRank_F;
        int mJumpingRank_B;

        // We need to mark whether nodes are atached or detached
        bool mAttached;
        
        // how many times have we seen this guy?
        int mCoverage;
        
        // is this a forward facing node?
        bool mIsForward;

        // we need to know which reads produced these nodes
        std::vector<StringToken> mReadHeaders;  // headers of all reads which contain these spacers
        ReadList mReadHolders;					// waste of the last var,  shut up.
};

#endif //CrisprNode_h
