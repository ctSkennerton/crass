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
#include "crassDefines.h"
#include "CrisprNode.h"
#include "StringCheck.h"

// we hash together string tokens to make a unique key for each spacer
typedef unsigned int SpacerKey;

inline SpacerKey makeSpacerKey(StringToken backST, StringToken frontST)
{
    //-----
    // make a spacer key from two string tokens
    //
	if(backST < frontST)
	{
		return (backST * 10000000) + frontST;
	}
	return (frontST * 10000000) + backST;
}

class SpacerInstance {
    public:
        SpacerInstance (void) {
            SI_SpacerSeqID = 0;
            SI_LeadingNode = NULL;
            SI_LastNode = NULL;
            SI_InstanceCount = 0;
            mSpacerRank = 0;
            mContigID = 0;
            mAttached = false;
        }
        
        SpacerInstance (StringToken spacerID);
        SpacerInstance (StringToken spacerID, CrisprNode * leadingNode, CrisprNode * lastNode);
        ~SpacerInstance () {}
        
        //
        // get / set
        //
        inline void incrementCount(void) { SI_InstanceCount++; }
        inline unsigned int getCount(void) { return SI_InstanceCount; }
        inline StringToken getID(void) { return SI_SpacerSeqID; }
        inline CrisprNode * getLeader(void) { return SI_LeadingNode; }
        inline CrisprNode * getLast(void) { return SI_LastNode; }
        inline bool isAttached(void) { return mAttached; }
        inline void setAttached(bool attached) { mAttached = attached; }
        
        //
        // contig functions
        //
        inline int getContigID(void) { return mContigID; }
        inline void setContigID(int CID) { mContigID = CID; }
        inline int getSpacerRank(void) { return mSpacerRank; }
        inline void setSpacerRank(int rank) { mSpacerRank = rank; }
        
    private:
        StringToken SI_SpacerSeqID;               // the StringToken of this spacer
        CrisprNode * SI_LeadingNode;              // the first node of this spacer
        CrisprNode * SI_LastNode;                 // the last node
        unsigned int SI_InstanceCount;            // the number of times this exact instance has been seen
        bool mAttached;							  // is this spacer attached?
        int mSpacerRank;						  // how many spacers come off this guy?
        int mContigID;							  // contig ID
};


#endif //SpacerInstance_h
