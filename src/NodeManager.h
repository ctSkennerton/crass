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
#include <map>
#include <string>
#include <fstream>
#include <queue>

// local includes
#include "NodeManager.h"
#include "crassDefines.h"
#include "CrisprNode.h"
#include "SpacerInstance.h"
#include "libcrispr.h"
#include "StringCheck.h"
#include "ReadHolder.h"
#include "GraphDrawingDefines.h"
#include "Rainbow.h"
#include "CrassXML.h"

// typedefs
typedef std::map<StringToken, CrisprNode *> NodeList;
typedef std::map<StringToken, CrisprNode *>::iterator NodeListIterator;

typedef std::map<SpacerKey, SpacerInstance *> SpacerList;
typedef std::map<SpacerKey, SpacerInstance *>::iterator SpacerListIterator;

typedef std::vector<CrisprNode *> NodeVector;
typedef std::vector<CrisprNode *>::iterator NodeVectorIterator;

typedef std::pair<CrisprNode *, CrisprNode *> CrisprNodePair;
typedef std::pair<CrisprNode *, CrisprNode *> CrisprNodePairIterator;

typedef std::vector<SpacerKey> SpacerVector;
typedef std::vector<SpacerKey>::iterator SpacerVectorIterator;

typedef std::map<int, SpacerVector *>ContigList;
typedef std::map<int, SpacerVector *>::iterator ContigListIterator;

//macros
#define makeKey(i,j) (i*100000)+j

class WalkingManager {
    CrisprNodePair WM_WalkingElem;
    bool WM_Direction;
    EDGE_TYPE WM_Wanted_Edge_Type;
    
	public:
	
		WalkingManager(CrisprNode * firstNode, CrisprNode * secondNode, EDGE_TYPE incommingEdge)
		{
			WM_WalkingElem.first = firstNode;
			WM_WalkingElem.second = secondNode;
			WM_Wanted_Edge_Type = incommingEdge;
		}
		WalkingManager(void){}
		~WalkingManager(void){}
	
	
		CrisprNode * getFirstNode(void) { return WM_WalkingElem.first; }
		CrisprNode * getSecondNode(void) { return WM_WalkingElem.second; }
		CrisprNodePair getWalkingElem(void) { return WM_WalkingElem; }
		EDGE_TYPE getEdgeType(void) { return WM_Wanted_Edge_Type; }
		void setFirstNode(CrisprNode * fn) { WM_WalkingElem.first = fn; }
		void setSecontNode(CrisprNode * sn) { WM_WalkingElem.second = sn; }
		void setWalkingElem(CrisprNodePair np) { WM_WalkingElem = np; }
		void setWantedEdge(EDGE_TYPE e) { WM_Wanted_Edge_Type = e; }
};

class NodeManager {
    public:

        NodeManager(std::string drSeq, const options * userOpts);
        ~NodeManager(void);

		bool addReadHolder(ReadHolder * RH);

        NodeListIterator nodeBegin(void) { return NM_Nodes.begin(); } 
        NodeListIterator nodeEnd(void) { return NM_Nodes.end(); }
        
    // get / set
    
        inline StringCheck * getStringCheck(void) { return &NM_StringCheck; }
		void findCapNodes(NodeVector * capNodes);                               // go through all the node and get a list of pointers to the nodes that have only one edge
		void findAllNodes(NodeVector * allNodes);
		void findAllNodes(NodeVector * capNodes, NodeVector * otherNodes);
		int findCapsAt(NodeVector * capNodes, bool searchForward, bool isInner, bool doStrict, CrisprNode * queryNode);
        EDGE_TYPE getOppositeEdgeType(EDGE_TYPE currentEdgeType);
        int getSpacerCount( bool showDetached);

    // Walking
		void walk(void);
		bool getEdge(WalkingManager * walkElem, CrisprNode * node, EDGE_TYPE * et);
		bool stepForType(WalkingManager * walkElem, EDGE_TYPE * et, CrisprNode ** detatchDelay);
		
    // Cleaning
        int cleanGraph(void);
        void clearBubbles(CrisprNode * rootNode, EDGE_TYPE currentEdgeType);
    
    // Contigs
        void findSpacerForContig(SpacerInstanceVector * sv, int contigID);
		int cleanSpacerGraph(void);
		void removeSpacerBubbles(void);
		int splitIntoContigs(void);
		void findAllForwardAttachedNodes(NodeVector * nodes);
		int buildSpacerGraph(void);
        void clearContigs(void);
        void contigiseForwardSpacers(std::queue<SpacerInstance *> * walkingQueue, SpacerInstance * SI);
        bool getForwardSpacer(SpacerInstance ** retSpacer, SpacerInstance * SI);
        bool getPrevSpacer(SpacerInstance ** retSpacer, SpacerInstance * SI);

    // Making purdy colours
        void setDebugColourLimits(void);
        void setSpacerColourLimits(void);

    // Printing / IO
        void printDebugGraph(std::ostream &dataOut, std::string title, bool showDetached, bool printBackEdges, bool longDesc);         // Print a graphviz style graph of the DRs and spacers
		void printDebugNodeAttributes(std::ostream& dataOut, CrisprNode * currCrisprNode, std::string colourCode, bool longDesc); 
        void printSpacerGraph(std::ostream &dataOut, std::string title, bool longDesc, bool showSingles);         // Print a graphviz style graph of the FULL spacers
        void printSpacerGraph(void);         // Print a graphviz style graph of the FULL spacers
        std::string getSpacerGraphLabel(SpacerInstance * spacer, bool longDesc);

    
        void printSpacerKey(std::ostream &dataOut, int numSteps, std::string groupNumber);         													// make a key for the spacer graph 	
        void dumpReads(std::string readsFileName, bool showDetached, bool split);												// dump reads to this file
        
    // XML
        void printXML(std::ofstream * XMLFile, int GID, bool showDetached);						// print this node managers portion of the XML file 
        void addSpacersToDOM(CrassXML * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached);
        void printAssemblyToDOM(CrassXML * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached);


    // Spacer dictionaries
        void printAllSpacers(void);
		void dumpSpacerDict(std::string spacerFileName, bool showDetached);

    private:
		
	// functions
		bool splitReadHolder(ReadHolder * RH);
        void addCrisprNodes(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void addSecondCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void addFirstCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void setUpperAndLowerCoverage(void);
     
    // members
        std::string NM_DirectRepeatSequence;  				// the sequence of this managers direct repeat
        NodeList NM_Nodes;                    				// list of CrisprNodes this manager manages
        SpacerList NM_Spacers;                				// list of all the spacers
        ReadList NM_ReadList;                 				// list of readholders
        StringCheck NM_StringCheck;           				// string check object for unique strings 
        Rainbow NM_DebugRainbow;              				// the Rainbow class for making colours
        Rainbow NM_SpacerRainbow;      				        // the Rainbow class for making colours
        const options * NM_Opts;              				// pointer to the user options structure
        int NM_NextContigID;									// next free contig ID (doubles as a counter)
        ContigList NM_Contigs; 								// our contigs
};


#endif // NodeManager_h
