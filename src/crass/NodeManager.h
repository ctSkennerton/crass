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
#include "Xml.h"
#include "StatsManager.h"


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
    SpacerInstancePair WM_WalkingElem;
    SI_EdgeDirection WM_Wanted_Edge_Type;
    
	public:
	
		WalkingManager(SpacerInstance * firstNode, SpacerInstance * secondNode, SI_EdgeDirection incommingEdge)
		{
			WM_WalkingElem.first = firstNode;
			WM_WalkingElem.second = secondNode;
			WM_Wanted_Edge_Type = incommingEdge;
		}
		WalkingManager(void){}
		~WalkingManager(void){}
	
	
		SpacerInstance * getFirstNode(void) { return WM_WalkingElem.first; }
		SpacerInstance * getSecondNode(void) { return WM_WalkingElem.second; }
        SpacerInstance * first(void) { return WM_WalkingElem.first; }
        SpacerInstance * second(void) { return WM_WalkingElem.second; }
    
		SpacerInstancePair getWalkingElem(void) { return WM_WalkingElem; }
		SI_EdgeDirection getEdgeType(void) { return WM_Wanted_Edge_Type; }
		void setFirstNode(SpacerInstance * fn) { WM_WalkingElem.first = fn; }
		void setSecondNode(SpacerInstance * sn) { WM_WalkingElem.second = sn; }
		void setWalkingElem(SpacerInstancePair np) { WM_WalkingElem = np; }
		void setWantedEdge(SI_EdgeDirection e) { WM_Wanted_Edge_Type = e; }
    
    SpacerInstance * shift(SpacerInstance * newNode);
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
        int getSpacerCountAndStats( bool showDetached);
        inline bool haveAnyFlankers(void){return (0 != NM_FlankerNodes.size());}

    // Walking
        bool getSpacerEdgeFromCap(WalkingManager * walkElem, SpacerInstance * nextSpacer);
        bool getSpacerEdgeFromCross(WalkingManager * walkElem, SpacerInstance * nextSpacer);
        bool stepThroughSpacerPath(WalkingManager * walkElem, SpacerInstance ** previousNode);
        bool walkFromCross(SpacerInstanceList * crossNodes);

    // Cleaning
        int cleanGraph(void);
        bool clearBubbles(CrisprNode * rootNode, EDGE_TYPE currentEdgeType);
    
    // Contigs
        void getAllSpacerCaps(SpacerInstanceVector * sv);
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
        bool printSpacerGraph(std::string& outFileName, std::string title, bool longDesc, bool showSingles);         // Print a graphviz style graph of the FULL spacers
        bool printSpacerGraph(void);         // Print a graphviz style graph of the FULL spacers
        std::string getSpacerGraphLabel(SpacerInstance * spacer, bool longDesc);

    
        void printSpacerKey(std::ostream &dataOut, int numSteps, std::string groupNumber);         													// make a key for the spacer graph 	
        void dumpReads(std::string readsFileName, bool showDetached, bool split);												// dump reads to this file
        
    // XML
        void printXML(std::ofstream * XMLFile, int GID, bool showDetached);						// print this node managers portion of the XML file 
        void addSpacersToDOM(crispr::XML * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached);
        void addFlankersToDOM(crispr::XML * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached);
        void printAssemblyToDOM(crispr::XML * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached);


    // Spacer dictionaries
        void printAllSpacers(void);
    
    // Flankers
        void generateFlankers(bool showDetached=false);

    private:
		
	// functions
		bool splitReadHolder(ReadHolder * RH);
        void addCrisprNodes(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void addSecondCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void addFirstCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt);
        void setContigIDForSpacers(SpacerInstanceVector * currentContigNodes);
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
        int NM_NextContigID;								// next free contig ID (doubles as a counter)
        ContigList NM_Contigs; 								// our contigs
        StatsManager<std::vector<size_t> > NM_SpacerLenStat;   // Keep a check on all of the spacer lengths for deciding whecher thay are a flanker or not
        SpacerInstanceVector NM_FlankerNodes;               // a list of spacers that are also flankers -- used only in the print functions
};


#endif // NodeManager_h
