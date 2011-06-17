// File: NodeManager.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Class to handle all the lists of nodes and class instances.
// Each cannonical form of direct repeat gets it's "own" NodeManager
// Thus a cannonical DR is a crispr is a NodeManager
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

#ifndef NodeManager_h
    #define NodeManager_h

// system includes
#include <iostream>
#include <vector>
#include <string>

// local includes
#include "NodeManager.h"
#include "crass_defines.h"
#include "CrisprNode.h"
#include "SpacerInstance.h"

// typedefs
typedef std::vector<CrisprNode *> NodeList;
typedef std::vector<CrisprNode *>::iterator NodeListIterator;
typedef std::vector<SpacerInstance *> SpacerList;
typedef std::vector<SpacerInstance *>::iterator SpacerListIterator;

class NodeManager {
    public:
        NodeManager(std::string drSeq);
        ~NodeManager(void);

    private:
        std::string mDirectRepeatSequence;  // the sequence of this managers direct repeat
        NodeList mNodes;                    // list of CrisprNodes this manager manages
        SpacerList mSpacers;                // list of all the spacers
};

#endif // NodeManager_h