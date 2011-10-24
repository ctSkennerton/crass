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
#include <map>
#include <string>
#include <sstream>
#include <fstream>

// local includes
#include <config.h>
#include "NodeManager.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "CrisprNode.h"
#include "SpacerInstance.h"
#include "libcrispr.h"
#include "StringCheck.h"
#include "ReadHolder.h"
#include "StringCheck.h"
#include "GraphDrawingDefines.h"
#include "Rainbow.h"

NodeManager::NodeManager(std::string drSeq, const options * userOpts)
{
    //-----
    // constructor
    //
    mDirectRepeatSequence = drSeq;
    mMinCoverage = 1000000;
    mMaxCoverage = 0;
    mOpts = userOpts;
    mStringCheck.setName("NM_" + drSeq);
}

NodeManager::~NodeManager(void)
{
    //-----
    // destructor
    //
    
    // clean up all the cripsr nodes
    NodeListIterator node_iter = mNodes.begin();
    while(node_iter != mNodes.end())
    {
        if(NULL != node_iter->second)
        {
            delete node_iter->second;
            node_iter->second = NULL;
        }
        node_iter++;
    }
    mNodes.clear();
    
    SpacerListIterator spacer_iter = mSpacers.begin();
    while(spacer_iter != mSpacers.end())
    {
        if(NULL != spacer_iter->second)
        {
            delete spacer_iter->second;
            spacer_iter->second = NULL;
        }
        spacer_iter++;
    }
    mSpacers.clear();

}

bool NodeManager::addReadHolder(ReadHolder * RH)
{
    //-----
    // add a readholder to this mofo
    //
    if (splitReadHolder(RH))
    {
        mReadList.push_back(RH);
        return true;
    }
    else
    {
        logError("Unable to split ReadHolder");
        return false;
    }
}

//----
// private function called from addReadHolder to split the read into spacers and pass it through to others
//
bool NodeManager::splitReadHolder(ReadHolder * RH)
{
    //-----
    // Split down a read holder and make some nodes
    //
    std::string working_str;
    CrisprNode * prev_node = NULL;
    
//    logInfo("Adding: " << RH->getHeader(),1);
    
    // add the header of this read to our stringcheck
    StringToken header_st = mStringCheck.addString(RH->getHeader());

    if(RH->getFirstSpacer(&working_str))
    {
        // do we have a direct repeat from the very beginning
        if (RH->startStopsAt(0) == 0) 
        {
            addCrisprNodes(&prev_node, working_str, header_st);
        } 
        else 
        {
            // we only want to add the second kmer, since it is anchored by the direct repeat
            addSecondCrisprNode(&prev_node, working_str, header_st);
        }
        
        // get all the spacers in the middle
        //check to see if we end with a direct repeat or a spacer
        if (RH->getSeqLength() == (int)RH->back() + 1) 
        {
            // direct repeat goes right to the end of the read take both
            while (RH->getNextSpacer(&working_str)) 
            {
                addCrisprNodes(&prev_node, working_str, header_st);
            }
        } 
        else 
        {
            // we end with an overhanging spacer so we want to break from the loop early
            // so that on the final time we only cut the first kmer            
            while (RH->getLastSpacerPos() < (int)RH->getStartStopListSize() - 1) 
            {
                //std::cout<<RH->getLastSpacerPos()<<" : "<<(int)RH->getStartStopListSize() - 1<<" : "<<working_str<<std::endl;
                RH->getNextSpacer(&working_str);
                addCrisprNodes(&prev_node, working_str, header_st);
            } 
            
            // get our last spacer
            if (RH->getNextSpacer(&working_str)) 
            {
                //std::cout<<working_str<<std::endl;
                addFirstCrisprNode(&prev_node, working_str, header_st);
            } 
        }
    }
    else
    {
        logError("Could not get a spacer for the read");
        return false;
    }
    return true;
}

//----
// Private function called from splitReadHolder to cut the kmers and make the nodes
//
void NodeManager::addCrisprNodes(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt)
{
    //-----
    // Given a spacer string, cut kmers from each end and make crispr nodes
    //
    // now cut kmers on either side and add the pair into the node manager 
    if (workingString.length() < CRASS_DEF_NODE_KMER_SIZE)
        return;
    
    std::string first_kmer = workingString.substr(0, CRASS_DEF_NODE_KMER_SIZE);
    std::string second_kmer = workingString.substr(workingString.length() - CRASS_DEF_NODE_KMER_SIZE, CRASS_DEF_NODE_KMER_SIZE );
    //std::cout<<"B: "<<first_kmer<<" : "<<second_kmer<<std::endl;
    
    CrisprNode * first_kmer_node;
    CrisprNode * second_kmer_node;
    
    // check to see if these kmers are already stored
    StringToken st1 = mStringCheck.getToken(first_kmer);
    StringToken st2 = mStringCheck.getToken(second_kmer);

    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = mStringCheck.addString(first_kmer);
        first_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        mNodes[st1] = first_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        first_kmer_node = mNodes[st1];
        (mNodes[st1])->incrementCount();
    }
    
    if(0 == st2)
    {
        st2 = mStringCheck.addString(second_kmer);

        second_kmer_node = new CrisprNode(st2);
        second_kmer_node->setForward(false);
        mNodes[st2] = second_kmer_node;
    }
    else
    {
        second_kmer_node = mNodes[st2];
        (mNodes[st2])->incrementCount();
    }
    
    // add in the read headers for the two CrisprNodes
    first_kmer_node->addReadHeader(headerSt);
    second_kmer_node->addReadHeader(headerSt);

    // the first kmers pair is the previous node which lay before it therefore bool is true
    // make sure prevNode is not NULL
    if (*prevNode != NULL/* && !seen_first*/) 
    {
        (*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
        first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
    }

    // now it's time to add the spacer
    SpacerInstance * curr_spacer;

    // check to see if we already have it here
    SpacerKey this_sp_key = makeSpacerKey(st1, st2);

    if(mSpacers.find(this_sp_key) == mSpacers.end())
    {
        // new instance
        StringToken sp_str_token = mStringCheck.addString(workingString);
        curr_spacer = new SpacerInstance(sp_str_token, first_kmer_node, second_kmer_node);
        mSpacers[this_sp_key] = curr_spacer;

        // make the inner edge
        first_kmer_node->addEdge(second_kmer_node, CN_EDGE_FORWARD);
        second_kmer_node->addEdge(first_kmer_node, CN_EDGE_BACKWARD);
    }
    else
    {
        // increment the number of times we've seen this guy
        (mSpacers[this_sp_key])->incrementCount();
    }

    *prevNode = second_kmer_node;
}

void NodeManager::addSecondCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt)
{
    if (workingString.length() < CRASS_DEF_NODE_KMER_SIZE)
        return;

    std::string second_kmer = workingString.substr(workingString.length() - CRASS_DEF_NODE_KMER_SIZE, CRASS_DEF_NODE_KMER_SIZE );
    CrisprNode * second_kmer_node;
    
    // check to see if these kmers are already stored
    StringToken st1 = mStringCheck.getToken(second_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = mStringCheck.addString(second_kmer);
        second_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        mNodes[st1] = second_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        second_kmer_node = mNodes[st1];
        (mNodes[st1])->incrementCount();
    }
    
    // add in the read headers for the this CrisprNode
    second_kmer_node->addReadHeader(headerSt);
    
    // add this guy in as the previous node for the next iteration
    *prevNode = second_kmer_node;

}

void NodeManager::addFirstCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt)
{
    if (workingString.length() < CRASS_DEF_NODE_KMER_SIZE)
        return;

    std::string first_kmer = workingString.substr(0, CRASS_DEF_NODE_KMER_SIZE);
    CrisprNode * first_kmer_node;
    //std::cout<<"F: "<<first_kmer<<std::endl;
    
    // check to see if these kmers are already stored
    StringToken st1 = mStringCheck.getToken(first_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = mStringCheck.addString(first_kmer);
        first_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        mNodes[st1] = first_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        first_kmer_node = mNodes[st1];
        (mNodes[st1])->incrementCount();
    }
    
    // add in the read headers for the this CrisprNode
    first_kmer_node->addReadHeader(headerSt);

    // check to see if we already have it here
    if(NULL != *prevNode)
    {
		SpacerKey this_sp_key = makeSpacerKey(st1, (*prevNode)->getID());
		if(mSpacers.find(this_sp_key) == mSpacers.end())
		{
			(*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
			first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
		}
    }
}

// Walking

void NodeManager::walk(void)
{
    // get the cap nodes
    NodeVector cap_nodes;
    findCapNodes(&cap_nodes);
    
    CrisprNode * detatch_delay = NULL;
    NodeVector detatch_list;
    
    WalkingManager * walk_elem = new WalkingManager();
    
    for (int max_step_cutoff = 1; max_step_cutoff <= CRASS_DEF_MAX_CLEANING; max_step_cutoff++) 
    {
        std::cout<<"min step cutoff: "<<max_step_cutoff<<std::endl;
        NodeVector::iterator cap_node_iter = cap_nodes.begin();

        while (cap_node_iter != cap_nodes.end()) 
        {
            int walk_steps = 0;
            EDGE_TYPE wanted_edge;
            if (getEdge(walk_elem, *cap_node_iter, &wanted_edge)) 
            {
                do {                    
                    if (detatch_delay != NULL) 
                    {
                        detatch_list.push_back(detatch_delay);
                    }
                    if (walk_steps > max_step_cutoff) 
                    {
                        break;
                    }
                    std::cout<<"Current node: "<<(walk_elem->getFirstNode())->getID()<<" -- "<< (walk_elem->getSecondNode())->getID()<<std::endl;
                    walk_steps++;

                } while (stepForType(walk_elem, &wanted_edge, &detatch_delay));
                std::cout<<"Nodes traversed in this walk: "<<walk_steps<<std::endl;
                
                if (walk_steps <= max_step_cutoff) 
                {
                    std::cout<<"Detatching nodes: ";
                    // detach nodes
                    NodeVector::iterator detach_vec_iter = detatch_list.begin();
                    while (detach_vec_iter != detatch_list.end()) 
                    {
                        std::cout<<(*detach_vec_iter)->getID()<<",";
                        (*detach_vec_iter)->detachNode();
                        detach_vec_iter++;
                    }
                    std::cout<<std::endl;
                } 
                // clear node vector
                detatch_list.clear();
            }
            cap_node_iter++;
        }
    }


}

bool NodeManager::getEdge(WalkingManager * walkElem, CrisprNode * currNode, EDGE_TYPE * et)
{
    // initaliser for the walking element 
    // from the cap node get the only edge and set that in the walking element
    
    
    // if there is a single forward jumping edge or single backward
    // jumping edge the node is a dead branch and needs to be cleaned
    if (((currNode)->getEdges(CN_EDGE_JUMPING_F))->size() == 1)
    {
        currNode->detachNode();
        return false;
    }
    
    else if((currNode->getEdges(CN_EDGE_JUMPING_B)->size()) == 1)
    {
        currNode->detachNode();
        return false;
    }
    
    // so if there is a single inner edge we've got our walking 
    // element but now we need to get the jumping edge equivelent
    else if((currNode->getEdges(CN_EDGE_FORWARD)->size()) == 1)
    {
        walkElem->setWantedEdge(CN_EDGE_JUMPING_F);
        *et = CN_EDGE_FORWARD;
    }

    else if((currNode->getEdges(CN_EDGE_BACKWARD)->size()) == 1)
    {
        walkElem->setWantedEdge(CN_EDGE_JUMPING_B);
        *et = CN_EDGE_BACKWARD;
    }
    else
    {
        logError("Trying to initialize a walking element using a node that is not a cap");
        return false;
    }    

    // get the second node for the walking element
    edgeListIterator edge_iter = (currNode->getEdges(*et))->begin();
    walkElem->setFirstNode(currNode);
    walkElem->setSecontNode(edge_iter->first);
    
    
    //
    return true;
}

bool NodeManager::stepForType(WalkingManager * walkElem, EDGE_TYPE * et, CrisprNode ** detatchDelay)
{
    // get all of the edges for our wanted edge type from the second node in the walking element
    edgeList * edge = (walkElem->getSecondNode())->getEdges(walkElem->getEdgeType());
    
    switch (edge->size()) {
        case 1:
            // there is only a single edge for the wanted type so we have a linear length
            break;
            
        default:
            // there is more than a single edge therefore we have a cross node
            break;
    }
    
    if (edge->size() == 1) 
    {
        edgeListIterator edge_iter = edge->begin();
        // if the node is attached
        if (edge_iter->second) 
        {
            // add the first node of the current walking element to the detach delay pointer
            *detatchDelay = walkElem->getFirstNode();
            
            // add the new node to the walking element
            walkElem->setFirstNode(edge_iter->first);
            // want opposite of what came in
            switch (*et) 
            {
                case CN_EDGE_BACKWARD:
                    walkElem->setWantedEdge(CN_EDGE_JUMPING_B);
                    *et = CN_EDGE_JUMPING_B;
                    return true;
                    break;
                 case CN_EDGE_FORWARD:
                    walkElem->setWantedEdge(CN_EDGE_JUMPING_F);
                    *et = CN_EDGE_JUMPING_F;
                    return true;
                    break;
                case CN_EDGE_JUMPING_B:
                    walkElem->setWantedEdge(CN_EDGE_BACKWARD);
                    *et = CN_EDGE_BACKWARD;
                    return true;
                    break;
                case CN_EDGE_JUMPING_F:
                    walkElem->setWantedEdge(CN_EDGE_FORWARD);
                    *et = CN_EDGE_FORWARD;
                    return true;
                    break;
                default:
                    logError("Could not set a new wanted edge for the walking element");
                    return false;
                    break;
            }
        } 
        else 
        {
            // pointing to an unattached node what do i do?
            return false;
        }
    } 
    else 
    {
        // cross node should break?
        return false;
    }
}

void NodeManager::findCapNodes(NodeVector * capNodes)
{    
    NodeListIterator all_node_iter = mNodes.begin();
    //std::cout<<"finding cap nodes:"<<std::endl;
    while (all_node_iter != mNodes.end()) 
    {
        int count = 0;
        
        // get the count of all attached nodes
        countEdgesForType(&count, all_node_iter->second, CN_EDGE_FORWARD);
        countEdgesForType(&count, all_node_iter->second, CN_EDGE_BACKWARD);
        countEdgesForType(&count, all_node_iter->second, CN_EDGE_JUMPING_B);
        countEdgesForType(&count, all_node_iter->second, CN_EDGE_JUMPING_F);
        
        if (count == 1) 
        {
            //std::cout<<(all_node_iter->second)->getID()<<",";
            capNodes->push_back(all_node_iter->second);
        }
        
        all_node_iter++;
    }
    //std::cout<<std::endl;
}

void NodeManager::countEdgesForType(int * count, CrisprNode * currNode, EDGE_TYPE edgeType)
{
    edgeList * curr_list = currNode->getEdges(edgeType);
    edgeListIterator edge_iter = curr_list->begin();
    while (edge_iter != curr_list->end()) 
    {
        if(edge_iter->second)
        {
            (*count)++;
        }
        edge_iter++;
    }
}

void NodeManager::findAllNodes(NodeVector * allNodes)
{
	//-----
	// make a nodevector of all of the nodes!
	//
	NodeListIterator all_node_iter = mNodes.begin();
	while (all_node_iter != mNodes.end()) 
	{
		allNodes->push_back(all_node_iter->second);
		all_node_iter++;
	}
}
// Cleaning


int NodeManager::cleanGraph(void)
{
    //-----
    // Clean all the bits off the graph mofo!
    //
	std::multimap<std::string, StringToken> seen_map;
	NodeVector nv;
	
	// get the cap nodes
	findAllNodes(&nv);
	
	// go through once and build the multimap
	NodeVectorIterator nv_iter = nv.begin();
	while(nv_iter != nv.end())
	{
		StringToken ST =  (*nv_iter)->getID();
		std::vector<std::string> read_headers = (*nv_iter)->getReadHeaders(&mStringCheck);
		std::vector<std::string>::iterator rh_iter = read_headers.begin();
		while(rh_iter != read_headers.end())
		{
			seen_map.insert(std::pair<std::string, StringToken>(*rh_iter, ST));
			rh_iter++;
		}
		nv_iter++;
	}
	// now print
	nv_iter = nv.begin();
	while(nv_iter != nv.end())
	{
		std::cout << (*nv_iter)->getID() << " : " << mStringCheck.getString((*nv_iter)->getID()) << std::endl;
		std::vector<std::string> read_headers = (*nv_iter)->getReadHeaders(&mStringCheck);
		std::vector<std::string>::iterator rh_iter = read_headers.begin();
		while(rh_iter != read_headers.end())
		{
			std::pair<std::multimap<std::string, StringToken>::iterator, std::multimap<std::string, StringToken>::iterator> m_str_iter = seen_map.equal_range(*rh_iter);
			std::multimap<std::string, StringToken>::iterator it2 = m_str_iter.first;
			while(it2 != m_str_iter.second)
			{
				std::cout << "  [" << (*it2).first << ", " << (*it2).second << "]" << std::endl;
				it2++;
			}
			rh_iter++;
		}
		nv_iter++;
		std::cout << "-------------" << std::endl;
	}
	std::cout << "==============" << std::endl;
	
	return 0;
}

// Spacer dictionaries
void NodeManager::dumpSpacerDict(std::string spacerFileName)
{
	//-----
	// Dump a spacer dictionary to file
	//
	std::ofstream spacer_file;
	spacer_file	.open(spacerFileName.c_str());
    if (spacer_file.good()) 
    {
    	SpacerListIterator spacer_iter = mSpacers.begin();
    	while(spacer_iter != mSpacers.end())
    	{
    		SpacerInstance * SI = spacer_iter->second;
    		std::string spacer = mStringCheck.getString(SI->getID());
        	spacer_file << spacer << " : " << SI->getCount() << " : " << (SI->getLeader())->getID() << " : " << (SI->getLast())->getID() << std::endl;
    		spacer_iter++;
    	}
    	spacer_file.close();
    }

}


// Making purdy colours
void NodeManager::setUpperAndLowerCoverage(void)
{
    // loop through all of the nodes and determine the upper and lower dounds for our graph
    NodeListIterator nl_iter = nodeBegin();
    while (nl_iter != nodeEnd()) 
    {
        int coverage = (nl_iter->second)->getCoverage();
        if (coverage > mMaxCoverage) 
        {
            mMaxCoverage = coverage;
        }
        else if(coverage < mMinCoverage)
        {
            mMinCoverage = coverage;
        }
        nl_iter++;
    }
//-DDEBUG#ifdef DEBUG
    logInfo("Max Node Coverage: "<<mMaxCoverage<<" Min Node Coverage: "<<mMinCoverage<<std::endl,5);
//-DDEBUG#endif
}


void NodeManager::setColourLimits(void)
{
    //-----
    // Make the colurs needed for printing the graphviz stuff
    //
    setUpperAndLowerCoverage();
    mRainbow.setType(mOpts->graphColourType);
    if (mOpts->coverageBins != -1) 
    {
        mRainbow.setLimits(mMinCoverage, mMaxCoverage, mOpts->coverageBins);
    } 
    else 
    {
    mRainbow.setLimits(mMinCoverage,mMaxCoverage);
    }
}

// Printing / IO


void NodeManager::printGraph(std::ostream &dataOut, std::string title, bool showDetached, bool printBackEdges)
{
    //-----
    // Print a graphviz style graph of the DRs and spacers
    //
    setColourLimits();
    
    gvGraphHeader(dataOut, title);
    NodeListIterator nl_iter = nodeBegin();
    // first loop to print out the nodes
    while (nl_iter != nodeEnd()) 
    {
        // check whether we should print
        if((nl_iter->second)->isAttached() | showDetached)
        {
            printNodeAttributes(dataOut, nl_iter->second ,mRainbow.getColour((nl_iter->second)->getCoverage()));
        }
        nl_iter++;
    }
    nl_iter = nodeBegin();
    
    // and go through again to print the edges
    while (nl_iter != nodeEnd()) 
    {
        // check whether we should print
        if((nl_iter->second)->isAttached() | showDetached)
        {
        	std::stringstream ss;
        	ss << (nl_iter->second)->getID() << "_" << mStringCheck.getString((nl_iter->second)->getID());
        	std::string label = ss.str();
            //std::cout<<(nl_iter->second)->getID()<<" : "<<mStringCheck.getString( (nl_iter->second)->getID() )<<std::endl;
            (nl_iter->second)->printEdges(dataOut, &mStringCheck, label, showDetached, printBackEdges );
        }
        nl_iter++;
    }
    gvGraphFooter(dataOut)
}

void NodeManager::printNodeAttributes(std::ostream& dataOut, CrisprNode * currCrisprNode, std::string colourCode)
{
    // print the node declaration
	std::stringstream ss;
	ss << currCrisprNode->getID() << "_" << mStringCheck.getString(currCrisprNode->getID());
	std::string label = ss.str();
    if(currCrisprNode->isForward())
    {
        gvNodeF(dataOut,label,colourCode);
    }
    else
    {
        gvNodeB(dataOut,label,colourCode);
    }
}



















