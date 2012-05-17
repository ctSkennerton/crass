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
#include <queue>

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
#include "StlExt.h"
#include "Exception.h"



SpacerInstance * WalkingManager::shift(SpacerInstance * newNode)
{
    SpacerInstance * old_node = WM_WalkingElem.first;
    WM_WalkingElem.first = WM_WalkingElem.second;
    WM_WalkingElem.second = newNode;
    return old_node;
}

NodeManager::NodeManager(std::string drSeq, const options * userOpts)
{
    //-----
    // constructor
    //
    NM_DirectRepeatSequence = drSeq;
    NM_Opts = userOpts;
    NM_StringCheck.setName("NM_" + drSeq);
    NM_NextContigID = 0;
    
}

NodeManager::~NodeManager(void)
{
    //-----
    // destructor
    //
    
    // clean up all the cripsr nodes
    NodeListIterator node_iter = NM_Nodes.begin();
    while(node_iter != NM_Nodes.end())
    {
        if(NULL != node_iter->second)
        {
            delete node_iter->second;
            node_iter->second = NULL;
        }
        node_iter++;
    }
    NM_Nodes.clear();
    
    SpacerListIterator spacer_iter = NM_Spacers.begin();
    while(spacer_iter != NM_Spacers.end())
    {
        if(NULL != spacer_iter->second)
        {
            delete spacer_iter->second;
            spacer_iter->second = NULL;
        }
        spacer_iter++;
    }
    NM_Spacers.clear();
    
    // delete contigs;
    clearContigs();
}

bool NodeManager::addReadHolder(ReadHolder * RH)
{
    //-----
    // add a readholder to this mofo
    //
    if (splitReadHolder(RH))
    {
        NM_ReadList.push_back(RH);
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
	
	// add the header of this read to our stringcheck

	StringToken header_st = NM_StringCheck.addString(RH->getHeader());
#ifdef SEARCH_SINGLETON
    SearchCheckerList::iterator debug_iter = debugger->find(RH->getHeader());
    if ( debug_iter != debugger->end()) {
        // an interesting read
        debug_iter->second.nmtoken(header_st);
    }
#endif
	//MI std::cout << std::endl << "----------------------------------\n" << RH->getHeader() << std::endl; 
	//MI std::cout << RH->splitApartSimple() << std::endl;
	
	if(RH->getFirstSpacer(&working_str))
	{
		//MI std::cout << "first SP: " << working_str << std::endl;
		try {
			// do we have a direct repeat from the very beginning
			if (RH->startStopsAt(0) == 0) 
			{
				//MI std::cout << "both" << std::endl;
				addCrisprNodes(&prev_node, working_str, header_st, RH);
			} 
			else 
			{
				//MI std::cout << "sec" << std::endl;
				// we only want to add the second kmer, since it is anchored by the direct repeat
				addSecondCrisprNode(&prev_node, working_str, header_st, RH);
			}
			
			// get all the spacers in the middle
			//check to see if we end with a direct repeat or a spacer
			if (RH->getSeqLength() == (int)RH->back() + 1) 
			{
				// direct repeat goes right to the end of the read take both
				//MI std::cout << "DR until end" << std::endl;
				while (RH->getNextSpacer(&working_str)) 
				{		
					//MI std::cout << "SP: " << working_str << std::endl;
					addCrisprNodes(&prev_node, working_str, header_st, RH);
				}
			} 
			else 
			{
				//MI std::cout << "SP at end" << std::endl;
				// we end with an overhanging spacer so we want to break from the loop early
				// so that on the final time we only cut the first kmer            
				while (RH->getLastSpacerPos() < (int)RH->getStartStopListSize() - 1) 
				{
					//std::cout<<RH->getLastSpacerPos()<<" : "<<(int)RH->getStartStopListSize() - 1<<" : "<<working_str<<std::endl;
					RH->getNextSpacer(&working_str);
					//MI std::cout << "SP: " << working_str << std::endl;
					addCrisprNodes(&prev_node, working_str, header_st, RH);
				} 
				
				// get our last spacer
				if (RH->getNextSpacer(&working_str)) 
				{
					//std::cout<<working_str<<std::endl;
					//MI std::cout << "last SP: " << working_str << std::endl;
					addFirstCrisprNode(&prev_node, working_str, header_st, RH);
				} 
			}
		} catch (crispr::substring_exception& e) {
			std::cerr<<e.what()<<std::endl;
			exit(99);
		} catch (...) {
			std::cerr<<"an unknown exception has occurred "<<__FILE__<<" : "<<__LINE__<<" : "<<__PRETTY_FUNCTION__<<std::endl; 
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
void NodeManager::addCrisprNodes(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt, ReadHolder * RH)
{
    //-----
    // Given a spacer string, cut kmers from each end and make crispr nodes
    //
    // now cut kmers on either side and add the pair into the node manager 
    if ((int)workingString.length() < NM_Opts->cNodeKmerLength)
        return;
    
    std::string first_kmer = workingString.substr(0, NM_Opts->cNodeKmerLength);
    std::string second_kmer = workingString.substr(workingString.length() - NM_Opts->cNodeKmerLength, NM_Opts->cNodeKmerLength );
    
    CrisprNode * first_kmer_node;
    CrisprNode * second_kmer_node;
    SpacerKey this_sp_key;
    
    // check to see if these kmers are already stored
    StringToken st1 = NM_StringCheck.getToken(first_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = NM_StringCheck.addString(first_kmer);
        first_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        NM_Nodes[st1] = first_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        first_kmer_node = NM_Nodes[st1];
        first_kmer_node->incrementCount();
    }
    
    StringToken st2 = NM_StringCheck.getToken(second_kmer);
    if(0 == st2)
    {
        st2 = NM_StringCheck.addString(second_kmer);
        second_kmer_node = new CrisprNode(st2);
        second_kmer_node->setForward(false);
        NM_Nodes[st2] = second_kmer_node;
    }
    else
    {
        second_kmer_node = NM_Nodes[st2];
        second_kmer_node->incrementCount();
    }
    // add in the read headers for the two CrisprNodes
    first_kmer_node->addReadHeader(headerSt);
    second_kmer_node->addReadHeader(headerSt);
    first_kmer_node->addReadHolder(RH);
    second_kmer_node->addReadHolder(RH);
    
    // the first kmers pair is the previous node which lay before it therefore bool is true
    // make sure prevNode is not NULL
    if (NULL != *prevNode) 
    {
        this_sp_key = makeSpacerKey(st1, (*prevNode)->getID());
        if(NM_Spacers.find(this_sp_key) == NM_Spacers.end())
        {
            (*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
            first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
        }
    }
    
#ifdef SEARCH_SINGLETON
    SearchCheckerList::iterator debug_iter = debugger->find(NM_StringCheck.getString(headerSt));
    if (debug_iter != debugger->end()) {
        // interesting read
        debug_iter->second.addNode(st1);       
        debug_iter->second.addNode(st2);
    }
#endif
    // now it's time to add the spacer
    SpacerInstance * curr_spacer;
    
    // check to see if we already have it here
    this_sp_key = makeSpacerKey(st1, st2);
    
    if(NM_Spacers.find(this_sp_key) == NM_Spacers.end())
    {
        // new instance
        StringToken sp_str_token = NM_StringCheck.getToken(workingString);
    	if(0 == sp_str_token)
    	{
            sp_str_token = NM_StringCheck.addString(workingString);
    	}
        curr_spacer = new SpacerInstance(sp_str_token, first_kmer_node, second_kmer_node);
        NM_Spacers[this_sp_key] = curr_spacer;
#ifdef SEARCH_SINGLETON
        if (debug_iter != debugger->end()) {
            debug_iter->second.addSpacer(workingString);
        }
#endif

        // make the inner edge
        first_kmer_node->addEdge(second_kmer_node, CN_EDGE_FORWARD);
        second_kmer_node->addEdge(first_kmer_node, CN_EDGE_BACKWARD);
    }
    else
    {
        // increment the number of times we've seen this guy
        (NM_Spacers[this_sp_key])->incrementCount();
    }
    
    *prevNode = second_kmer_node;
}

void NodeManager::addSecondCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt, ReadHolder * RH)
{
    if ((int)workingString.length() < NM_Opts->cNodeKmerLength)
        return;
    
    std::string second_kmer = workingString.substr(workingString.length() - NM_Opts->cNodeKmerLength, NM_Opts->cNodeKmerLength );
    CrisprNode * second_kmer_node;
    
    // check to see if these kmers are already stored
    StringToken st2 = NM_StringCheck.getToken(second_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st2)
    {
        // first time we've seen this guy. Make some new objects
        st2 = NM_StringCheck.addString(second_kmer);
        second_kmer_node = new CrisprNode(st2);
        second_kmer_node->setForward(false);
        
        // add them to the pile
        NM_Nodes[st2] = second_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        second_kmer_node = NM_Nodes[st2];
        (NM_Nodes[st2])->incrementCount();
    }
#ifdef SEARCH_SINGLETON
    SearchCheckerList::iterator debug_iter = debugger->find(NM_StringCheck.getString(headerSt));
    if (debug_iter != debugger->end()) {
        // interesting read
        debug_iter->second.addNode(st2);
    }
#endif
    // add in the read headers for the this CrisprNode
    second_kmer_node->addReadHeader(headerSt);
    second_kmer_node->addReadHolder(RH);
    
    // add this guy in as the previous node for the next iteration
    *prevNode = second_kmer_node;
    
    // there is no one yet to make an edge
}

void NodeManager::addFirstCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt, ReadHolder * RH)
{
    if ((int)workingString.length() < NM_Opts->cNodeKmerLength)
        return;
    
    std::string first_kmer = workingString.substr(0, NM_Opts->cNodeKmerLength);
    CrisprNode * first_kmer_node;
    
    // check to see if these kmers are already stored
    StringToken st1 = NM_StringCheck.getToken(first_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = NM_StringCheck.addString(first_kmer);
        first_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        NM_Nodes[st1] = first_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        first_kmer_node = NM_Nodes[st1];
        (NM_Nodes[st1])->incrementCount();
    }
#ifdef SEARCH_SINGLETON
    SearchCheckerList::iterator debug_iter = debugger->find(NM_StringCheck.getString(headerSt));
    if (debug_iter != debugger->end()) {
        // interesting read
        debug_iter->second.addNode(st1);       
    }
#endif
    // add in the read headers for the this CrisprNode
    first_kmer_node->addReadHeader(headerSt);
    first_kmer_node->addReadHolder(RH);
    
    // check to see if we already have it here
    if(NULL != *prevNode)
    {
        SpacerKey this_sp_key = makeSpacerKey(st1, (*prevNode)->getID());
        if(NM_Spacers.find(this_sp_key) == NM_Spacers.end())
        {
            (*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
            first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
        }
    }
}

// Walking



void NodeManager::findCapNodes(NodeVector * capNodes)
{
    //-----
    // find cap nodes!
    //
    capNodes->clear();
    
    NodeListIterator all_node_iter = NM_Nodes.begin();
    while (all_node_iter != NM_Nodes.end()) 
    {
        if((all_node_iter->second)->isAttached())
        {
            if ((all_node_iter->second)->getTotalRank() == 1) 
            {
                capNodes->push_back(all_node_iter->second);
            }
        }
        all_node_iter++;
    }
}

void NodeManager::findAllNodes(NodeVector * allNodes)
{
    //-----
    // make a nodevector of all of the nodes!
    //
    allNodes->clear();
    
    NodeListIterator all_node_iter = NM_Nodes.begin();
    while (all_node_iter != NM_Nodes.end()) 
    {
        if((all_node_iter->second)->isAttached())
        {
            allNodes->push_back(all_node_iter->second);
            all_node_iter++;
        }
    }
}

void NodeManager::findAllNodes(NodeVector * capNodes, NodeVector * otherNodes)
{
    //-----
    // make nodevectors of all of the nodes, split between cap and other 
    //
    capNodes->clear();
    otherNodes->clear();
    
    NodeListIterator all_node_iter = NM_Nodes.begin();
    while (all_node_iter != NM_Nodes.end()) 
    {
        if((all_node_iter->second)->isAttached())
        {
            int rank = (all_node_iter->second)->getTotalRank(); 
            if (rank == 1) 
            { capNodes->push_back(all_node_iter->second); }
            else
            { otherNodes->push_back(all_node_iter->second); }
        }
        all_node_iter++;
    }
}

int NodeManager::findCapsAt(NodeVector * capNodes, bool searchForward, bool isInner, bool doStrict, CrisprNode * queryNode)
{
    //-----
    // Populate the passed nodevector with the capnodes at queryNode
    // and return the size of the vector
    //
    // if doStrict is true then the function will only return a non-zero result 
    // when the the query node is joined ONLY to cap nodes (of the given edge type)
    //
    capNodes->clear();
    
    if(queryNode->isAttached())
    {
        // first, find what type of edge we are searching for.
        edgeListIterator el_iter;
        edgeList * el;
        if(searchForward)
        {
            if(isInner)
                el = queryNode->getEdges(CN_EDGE_FORWARD);
            else
                el = queryNode->getEdges(CN_EDGE_JUMPING_F);
        }
        else
        {
            if(isInner)
                el = queryNode->getEdges(CN_EDGE_BACKWARD);
            else
                el = queryNode->getEdges(CN_EDGE_JUMPING_B);
        }
        el_iter = el->begin();
        while(el_iter != el->end())
        {
            // make sure we only look at attached edges
            if(el_iter->second)
            {
                CrisprNode * attached_node = el_iter->first;
                if(1 == attached_node->getTotalRank())
                {
                    // this guy is a cap!
                    capNodes->push_back(attached_node);
                }
                else
                {
                    // non-cap of the same edge type
                    if(doStrict)
                    {
                        capNodes->clear();
                        return 0;
                    }
                }
            }
            
            el_iter++;
        }
    }
    return (int)capNodes->size();
}


bool NodeManager::getSpacerEdgeFromCap(WalkingManager * walkElem, SpacerInstance * currentSpacer)
{
    if (currentSpacer->getSpacerRank() == 1) 
    {
        SpacerEdgeVector_Iterator iter = currentSpacer->begin();
        
        while (iter != currentSpacer->end()) 
        {
            if (((*iter)->edge)->isAttached()) 
            {
                if (0 == ((*iter)->edge)->getContigID()) 
                {
                    walkElem->setSecondNode((*iter)->edge);
                    walkElem->setFirstNode(currentSpacer);
                    walkElem->setWantedEdge((*iter)->d);
                } 
                else 
                {
                    currentSpacer->setContigID(((*iter)->edge)->getContigID());
                    return false;
                }
            } 
            else 
            {
                return false;
            }
            iter++;
        }
    }
    else 
    {
        return false;
    }
    if (walkElem->first() == NULL || walkElem->second() == NULL) {
        return false;
    } else {
        return true;
    }
}

bool NodeManager::getSpacerEdgeFromCross(WalkingManager * walkElem, SpacerInstance * currentSpacer )
{
    // check that the edge is a path node
    if (currentSpacer->getSpacerRank() == 2) 
    {
        SpacerEdgeVector_Iterator iter = currentSpacer->begin();
        
        while (iter != currentSpacer->end()) 
        {
            if (((*iter)->edge)->isAttached()) 
            {
                if (0 == ((*iter)->edge)->getContigID()) 
                {
                    walkElem->setSecondNode((*iter)->edge);
                    walkElem->setFirstNode(currentSpacer);
                    walkElem->setWantedEdge((*iter)->d);
                    return true;
                } 
            } 
            else 
            {
                return false;
            }
            iter++;
        }
    }
    else 
    {
        return false;
    }
    if (walkElem->first() == NULL || walkElem->second() == NULL) {
        return false;
    } else {
        return true;
    }
}

bool NodeManager::stepThroughSpacerPath(WalkingManager * walkElem, SpacerInstance ** previousNode)
{
    switch ((walkElem->second())->getSpacerRank()) 
    {
        case 2:
        {
            // path node?
            // check to see if the other edge is going in the same direction as wantedDirection
            SpacerEdgeVector_Iterator iter = (walkElem->second())->begin();
            while (iter != (walkElem->second())->end()) 
            {
                // make sure the edge is attached, in the right direction, not the incomming edge and not assigned to a contig already
                if ((*iter)->edge->isAttached()  && 
                    ((*iter)->d == walkElem->getEdgeType()) && 
                    ((*iter)->edge->getID() != walkElem->first()->getID()) && 
                    ((*iter)->edge->getContigID() == 0)) 
                {
                    *previousNode = walkElem->shift((*iter)->edge);
                    return true;
                } 
                iter++;
            }
            break;
        }
        case 1:
            //
        default:
        {
            return false;
            // cross node
            break;
        }
    }
    return false;
}


// Cleaning
int NodeManager::cleanGraph(void)
{
    //-----
    // Clean all the bits off the graph mofo!
    //
    // keep going while we're detaching stuff
    bool some_detached = true;
    
    while(some_detached)
    {
        std::multimap<CrisprNode *, CrisprNode *> fork_choice_map;
        NodeVector nv_cap, nv_other, detach_list;
        NodeVectorIterator nv_iter;
        some_detached = false;
        
        // get all the nodes
        findAllNodes(&nv_cap, &nv_other);
    
        // First do caps
        nv_iter = nv_cap.begin();
        while(nv_iter != nv_cap.end())
        {
            // we can just lop off caps joined by jumpers (perhaps)
            if ((*nv_iter)->getInnerRank() == 0)
            {
                // make sure that this guy is linked to a cross node
                edgeList * el;
                if(0 != (*nv_iter)->getRank(CN_EDGE_JUMPING_F))
                    el = (*nv_iter)->getEdges(CN_EDGE_JUMPING_F);
                else
                    el = (*nv_iter)->getEdges(CN_EDGE_JUMPING_B);
                
                // there is only one guy in this list!
                int other_rank = ((el->begin())->first)->getTotalRank();
                if(other_rank != 2)
                    detach_list.push_back(*nv_iter);
            }
            else
            {
                // make sure that this guy is linked to a cross node
                edgeList * el;
                bool is_forward;
                if(0 != (*nv_iter)->getRank(CN_EDGE_FORWARD))
                {
                    el = (*nv_iter)->getEdges(CN_EDGE_FORWARD);
                    is_forward = false;
                }
                else
                {
                    el = (*nv_iter)->getEdges(CN_EDGE_BACKWARD);
                    is_forward = true;
                }
                
                // there is only one guy in this list!
                CrisprNode * joining_node = ((el->begin())->first); 
                int other_rank = joining_node->getTotalRank();
                if(other_rank != 2)
                {
                    // this guy joins onto a crossnode
                    // check to see if he is the only cap here!
                    NodeVector caps_at_join;
                    if(findCapsAt(&caps_at_join, is_forward, true, true, joining_node) > 1)
                    {
                        // this is a fork at the end of an arm
                        fork_choice_map.insert(std::pair<CrisprNode *, CrisprNode *>(joining_node, *nv_iter)) ;
                    }
                    else
                    {
                        // the only cap at a cross. NUKE!
                        detach_list.push_back(*nv_iter);
                    }
                }
            }
            nv_iter++;
        }
        
        // make coverage decisions for end forks
        std::map<CrisprNode *, int> best_coverage_map_cov;
        std::map<CrisprNode *, CrisprNode *> best_coverage_map_node;
        std::multimap<CrisprNode *, CrisprNode *>::iterator fcm_iter = fork_choice_map.begin();
        while(fcm_iter != fork_choice_map.end())
        {
            if(best_coverage_map_cov.find(fcm_iter->first) == best_coverage_map_cov.end())
            {
                // first one!
                best_coverage_map_cov.insert(std::pair<CrisprNode *, int>(fcm_iter->first, ((*fcm_iter).second)->getCoverage()));
                best_coverage_map_node.insert(std::pair<CrisprNode *, CrisprNode *>(fcm_iter->first, (*fcm_iter).second));
            }
            else
            {
                // one is already here!
                int new_cov = ((*fcm_iter).second)->getCoverage();
                if(best_coverage_map_cov[fcm_iter->first] < new_cov)
                {
                    // the new one is better!
                    best_coverage_map_cov[fcm_iter->first] = ((*fcm_iter).second)->getCoverage();
                    best_coverage_map_node[fcm_iter->first] = (*fcm_iter).second;
                }
            }
            fcm_iter++;
        }
        fcm_iter = fork_choice_map.begin();
        while(fcm_iter != fork_choice_map.end())
        {
            if(best_coverage_map_node[fcm_iter->first] != (*fcm_iter).second)
            {
                // not the best one!
                detach_list.push_back((*fcm_iter).second);
            }
            fcm_iter++;
        }

        // check to see if we'll need to do this again
        if(detach_list.size() > 0)
            some_detached = true;
        
        // finally, detach!
        nv_iter = detach_list.begin();
        while(nv_iter != detach_list.end())
        {
            (*nv_iter)->detachNode();
            nv_iter++;
        }
    
        // refresh the node lists
        findAllNodes(&nv_cap, &nv_other);
    
        // then do bubbles
        nv_iter = nv_other.begin();
        while(nv_iter != nv_other.end())
        {
            switch ((*nv_iter)->getTotalRank()) 
            {
                case 2:
                {
                    // check that there is one inner and one jumping edge
                    if (!((*nv_iter)->getInnerRank() && (*nv_iter)->getJumpingRank())) 
                    {
    #ifdef DEBUG
                        logInfo("node "<<(*nv_iter)->getID()<<" has only two edges of the same type -- cannot be linear -- detaching", 8);
    #endif
                        (*nv_iter)->detachNode();
                        some_detached = true;
                    }
                    break;
                }
                case 1:
                case 0:
                    break;
                default:
                {
                    // get the rank for the the inner and jumping edges.
                    if((*nv_iter)->getInnerRank() != 1)
                    {
                        // there are multiple inner edges for this guy
                        if(clearBubbles(*nv_iter, CN_EDGE_FORWARD))
                        	some_detached = true;
                    }
                    
                    if((*nv_iter)->getJumpingRank() != 1)
                    {
                        // there are multiple jumping edges for this guy
                        if(clearBubbles(*nv_iter, CN_EDGE_JUMPING_F))
                        	some_detached = true;
                    }
                    break;
                }
            }        
            nv_iter++;
        }
    }
    return 0;
}

bool NodeManager::clearBubbles(CrisprNode * rootNode, EDGE_TYPE currentEdgeType)
{
	//-----
	// Return true if something got detached
	//
	bool some_detached = false;
	
    // get a list of edges
    edgeList * curr_edges = rootNode->getEdges(currentEdgeType);
    
    // the key is the hashed values of both the root node and the edge
    // the value is the node id of the edge
    std::map<int, int> bubble_map;
    
    // now go through each of the edges and make a hashed key for the edge 
    edgeListIterator curr_edges_iter = curr_edges->begin();
    while(curr_edges_iter != curr_edges->end())
    {
        if ((curr_edges_iter->first)->isAttached()) 
        {
            // we want to go through all the edges of the nodes above (2nd degree separation)
            // and since we used the forward edges to get here we now want the opposite (Jummping_F)
            edgeList * edges_of_curr_edge = (curr_edges_iter->first)->getEdges(getOppositeEdgeType(currentEdgeType));
            
            edgeListIterator edges_of_curr_edge_iter = edges_of_curr_edge->begin();
            while (edges_of_curr_edge_iter != edges_of_curr_edge->end()) 
            {
                // make sue that this guy is attached
                if ((edges_of_curr_edge_iter->first)->isAttached()) 
                {
                    // so now we're at the second degree of separation for our edges
                    // again make a key but check to see if the key exists in the hash
                    
                    int new_key = makeKey(rootNode->getID(), (edges_of_curr_edge_iter->first)->getID());
                    if (bubble_map.find(new_key) == bubble_map.end()) 
                    {
                        // first time we've seen him
                        bubble_map[new_key] = (curr_edges_iter->first)->getID();
                    } 
                    else 
                    {
                        // aha! he is pointing back onto the same guy as someone else.  We have a bubble!
                        //get the CrisprNode of the first guy
                        
                        CrisprNode * first_node = NM_Nodes[bubble_map[new_key]];
#ifdef DEBUG
                        logInfo("Bubble found conecting "<<rootNode->getID()<<" : "<<first_node->getID()<<" : "<<(edges_of_curr_edge_iter->first)->getID(), 8);
#endif
                        //perform a coverage test on the nodes that end up here and kill the one with the least coverage
                        
                        // TODO: this is a pretty dumb test, since the coverage between the two nodes could be very similar
                        // for example one has a coverage of 20 and the other a coverage of 18.  The node with 18 will get
                        // removed but should it have been?  The way to fix this would be to create some infastructure in
                        // NodeManager to calculate the average and stdev of the coverage and then remove a node only if
                        // it is below 1 stdev of the average, else it could be a biological thing that this bubble exists.
                        
                        if (first_node->getDiscountedCoverage() > (curr_edges_iter->first)->getDiscountedCoverage()) 
                        {
                            // the first guy has greater coverage so detach our current node
                            (curr_edges_iter->first)->detachNode();
                            some_detached = true;
#ifdef DEBUG
                            logInfo("Detaching "<<(curr_edges_iter->first)->getID()<<" as it has lower coverage", 8);
#endif
                        } 
                        else 
                        {
                            // the first guy was lower so kill him
                            first_node->detachNode();
                            some_detached = true;
#ifdef DEBUG
                            logInfo("Detaching "<<first_node->getID()<<" as it has lower coverage", 8);
#endif
                            // replace the existing key (to check for triple bubbles)
                            bubble_map[new_key] = (curr_edges_iter->first)->getID();
                        }
                    }
                }
                edges_of_curr_edge_iter++;
            }
        }
        curr_edges_iter++;
    }
    return some_detached;
}

EDGE_TYPE NodeManager::getOppositeEdgeType(EDGE_TYPE currentEdgeType)
{
    switch (currentEdgeType) {
        case CN_EDGE_BACKWARD:
            return CN_EDGE_JUMPING_B;
            break;
        case CN_EDGE_FORWARD:
            return CN_EDGE_JUMPING_F;
            break;
        case CN_EDGE_JUMPING_B:
            return CN_EDGE_BACKWARD;
            break;
        case CN_EDGE_JUMPING_F:
            return CN_EDGE_FORWARD;
            break;
        default:
            logError("Cannot find the opposite edge type for input");
            return CN_EDGE_ERROR;
            break;
    }
}

int NodeManager::getSpacerCountAndStats(bool showDetached, bool excludeFlankers)
{
    int number_of_spacers = 0;
    SpacerListIterator sp_iter; NM_Spacers.begin();
    for (sp_iter = NM_Spacers.begin(); sp_iter != NM_Spacers.end(); ++sp_iter) 
    {
        SpacerInstance * current_spacer = sp_iter->second;
        if (showDetached || current_spacer->isAttached()) 
        {
            if (excludeFlankers & current_spacer->isFlanker()) {
                continue;
            }
            // add in some stats for the spacers
            std::string spacer = NM_StringCheck.getString((sp_iter->second)->getID());
            NM_SpacerLenStat.add(spacer.length());
            number_of_spacers++;
        }
    }
    return number_of_spacers;
}

// Contigs    
void NodeManager::clearContigs(void)
{
    //-----
    // Clear all contig information
    //
    ContigListIterator cl_iter = NM_Contigs.begin();
    while(cl_iter != NM_Contigs.end())
    {
        if(NULL != cl_iter->second)
        {
            SpacerVectorIterator sp_iter = (cl_iter->second)->begin();
            while(sp_iter != (cl_iter->second)->end())
            {
                (NM_Spacers[*sp_iter])->setContigID(0);
                sp_iter++;
            }
            delete cl_iter->second;
        }
        cl_iter++;
    }
    NM_NextContigID = 0;
}

void NodeManager::findAllForwardAttachedNodes(NodeVector * nodes)
{
    //-----
    // make nodevectors of all of the nodes, split between cap, cross, path and other 
    //
    nodes->clear();
    
    NodeListIterator all_node_iter = NM_Nodes.begin();
    while (all_node_iter != NM_Nodes.end()) 
    {
        if((all_node_iter->second)->isAttached()&& (all_node_iter->second)->isForward())
        {
            nodes->push_back(all_node_iter->second);
        }
        all_node_iter++;
    }
}

int NodeManager::buildSpacerGraph(void)
{
    //-----
    // For all forward nodes, count the number of ongoing spacers
    // make spacer edges if told to do so
    //
    SpacerListIterator spacers_iter = NM_Spacers.begin();
    while(spacers_iter != NM_Spacers.end()) 
    {
        // get the last node of this spacer
        CrisprNode * rq_leader_node = (spacers_iter->second)->getLeader();
        CrisprNode * rq_last_node = (spacers_iter->second)->getLast();
        
        if(rq_last_node->isAttached() && rq_leader_node->isAttached())
        {
            // mark this guy as attached
            (spacers_iter->second)->setAttached(true);
            
            // now get all the jumping forward edges from this node
            edgeList * qel = rq_last_node->getEdges(CN_EDGE_JUMPING_F);
            edgeListIterator qel_iter = qel->begin();
            while(qel_iter != qel->end())
            {
                if((qel_iter->first)->isAttached() && (qel_iter->first)->isForward())
                {
                    // a forward attached node. Now check for inner edges.
                    edgeList * el = (qel_iter->first)->getEdges(CN_EDGE_FORWARD);
                    edgeListIterator el_iter = el->begin();
                    while(el_iter != el->end())
                    {
                        if((el_iter->first)->isAttached())
                        {
                            // bingo!
                            SpacerInstance * next_spacer = NM_Spacers[makeSpacerKey((el_iter->first)->getID(), (qel_iter->first)->getID())];
                            
                            if (next_spacer == spacers_iter->second) {
                                //logError("Spacer "<<spacers_iter->second << " with id "<< (spacers_iter->second)->getID()<< " has an edge to itself... aborting edge "<<next_spacer <<" : "<< spacers_iter->second);
                            } 
                            else 
                            {
                                // we can add an edge for these two spacers
                                // add the forward edge to the next spacer
                                spacerEdgeStruct * new_edge = new spacerEdgeStruct();
                                new_edge->edge = next_spacer;
                                new_edge->d = FORWARD;
                                (spacers_iter->second)->addEdge(new_edge);
                                
                                // add the corresponding reverse edge to the current spacer
                                spacerEdgeStruct * new_edge2 = new spacerEdgeStruct();
                                new_edge2->edge = spacers_iter->second;
                                new_edge2->d = REVERSE;
                                next_spacer->addEdge(new_edge2);
                            }

                        }
                        el_iter++;
                    }
                }
                qel_iter++;
            }
        }
        else
        {
            (spacers_iter->second)->setAttached(false);    
        }
        spacers_iter++;
    }
    return 0;
}

void NodeManager::getAllSpacerCaps(SpacerInstanceVector * sv)
{
    SpacerListIterator sp_iter = NM_Spacers.begin();
    while (sp_iter != NM_Spacers.end()) 
    {
        if((sp_iter->second)->isAttached() )
        {
            if ((sp_iter->second)->getSpacerRank() == 1) 
            {
                sv->push_back(sp_iter->second);
            }
        }
        sp_iter++;
    }
}


void NodeManager::findSpacerForContig(SpacerInstanceVector * sv, int contigID)
{
    SpacerListIterator sp_iter = NM_Spacers.begin();
    
    // first we build up the contents of the walking queue
    while(sp_iter != NM_Spacers.end())
    {
        if((sp_iter->second)->isAttached() )
        {
            if ((sp_iter->second)->getContigID() == contigID) 
            {
                sv->push_back(sp_iter->second);
            }
        }
        sp_iter++;
    }
}

int NodeManager::cleanSpacerGraph(void)
{
    //-----
    // Clean up the spacer graph
    //
    int round  = 0;
    bool cleaned_some = true;
    while(cleaned_some)
    {
        round++;
        logInfo("Cleaning round: " << round, 2);
        cleaned_some = false;
        
        // remove fur
        SpacerListIterator sp_iter = NM_Spacers.begin();
        while(sp_iter != NM_Spacers.end())
        {
            if((sp_iter->second)->isAttached())
            {
                if(sp_iter->second->isFur())
                {
                    //std::cout << "a: " << sp_iter->second << std::endl;
                    sp_iter->second->detachFromSpacerGraph();
                    cleaned_some = true;
                }
            }
            sp_iter++;
        }
        
        // remove non-viable nodes
        sp_iter = NM_Spacers.begin();
        while(sp_iter != NM_Spacers.end())
        {
            if((sp_iter->second)->isAttached())
            {
                if(!sp_iter->second->isViable())
                {
                    //std::cout << "b: " << sp_iter->second << std::endl;
                    sp_iter->second->detachFromSpacerGraph();
                    cleaned_some = true;
                }
            }
            sp_iter++;
        }
        
        // remove bubbles
        removeSpacerBubbles();
    }
    return 0;
}

void NodeManager::removeSpacerBubbles(void)
{
    //-----
    // remove bubbles from the spacer graph
    //
    std::map<SpacerKey, SpacerInstance *> bubble_map;
    SpacerInstanceVector detach_list;
    SpacerListIterator sp_iter = NM_Spacers.begin();
    while(sp_iter != NM_Spacers.end())
    {
        if((sp_iter->second)->isAttached())
        {
            // we only car about rank 2 or over nodes
            if(2 <= (sp_iter->second)->getSpacerRank())
            {
                // first make a list of the forward and backward spacers
                SpacerEdgeVector_Iterator edge_iter = (sp_iter->second)->begin();
                SpacerInstanceVector f_spacers, r_spacers;
                while(edge_iter != (sp_iter->second)->end())
                {
                    if((*edge_iter)->d == REVERSE)
                        r_spacers.push_back((*edge_iter)->edge);
                    else
                        f_spacers.push_back((*edge_iter)->edge);
                    edge_iter++;
                }
                
                // now make a list of spacer keys 
                SpacerInstanceVector_Iterator r_edge_iter = r_spacers.begin();
                while(r_edge_iter != r_spacers.end())
                {
                    SpacerInstanceVector_Iterator f_edge_iter =  f_spacers.begin();
                    while(f_edge_iter != f_spacers.end())
                    {
                        // make a key
                        SpacerKey tmp_key = makeSpacerKey((*r_edge_iter)->getID(), (*f_edge_iter)->getID());
                        // check if we've seen this key before
                        std::map<SpacerKey, SpacerInstance *>::iterator bm_iter = bubble_map.find(tmp_key);
                        if(bm_iter == bubble_map.end())
                        {
                            // first time
                            bubble_map[tmp_key] = sp_iter->second;
                        }
                        else
                        {
                            // bubble! -- check the coverages!
                            if(bubble_map[tmp_key]->getCount() < (sp_iter->second)->getCount())
                            {
                                // stored guy has lower coverage!
                                detach_list.push_back(bubble_map[tmp_key]);
                                bubble_map[tmp_key] = sp_iter->second;
                            }
                            else if((sp_iter->second)->getCount() < bubble_map[tmp_key]->getCount())
                            {
                                // new guy has lower coverage!
                                detach_list.push_back(sp_iter->second);
                            }
                            else
                            {
                                // coverages are equal, kill the one with the lower rank
                                if(bubble_map[tmp_key]->getSpacerRank() < (sp_iter->second)->getSpacerRank())
                                {
                                    // stored guy has lower coverage!
                                    detach_list.push_back(bubble_map[tmp_key]);
                                    bubble_map[tmp_key] = sp_iter->second;
                                }
                                else
                                {
                                    // new guy has lower or equal coverage!
                                    detach_list.push_back(sp_iter->second);
                                }
                            }
                        }
                        f_edge_iter++;
                    }
                    r_edge_iter++;
                }
            }
        }
        sp_iter++;
    }
            
    // detach all on the detach list!
    SpacerInstanceVector_Iterator dl_iter = detach_list.begin();
    while(dl_iter != detach_list.end())
    {
        (*dl_iter)->detachFromSpacerGraph();
        dl_iter++;
    }
    
}

int NodeManager::splitIntoContigs(void)
{
    //-----
    // split the group into contigs 
    //
    
    
    // get all of the cap nodes
    SpacerInstanceVector start_walk_nodes;
    getAllSpacerCaps(&start_walk_nodes);
    
    SpacerInstanceList cross_nodes;
    
    WalkingManager * walk_elem = new WalkingManager();
    // walk from the cap nodes to a cross node
    SpacerInstanceVector_Iterator cap_node_iter = start_walk_nodes.begin();
    while (cap_node_iter != start_walk_nodes.end())
    {
        SpacerInstanceVector current_contig_spacers;
        NM_NextContigID++;
        if (getSpacerEdgeFromCap(walk_elem,*cap_node_iter))
        {
            SpacerInstance * previous_spacer = NULL;
            
            //current_contig_spacers.push_back(*cap_node_iter);
            do { 
                if (NULL != previous_spacer) 
                {
                    
                    current_contig_spacers.push_back(previous_spacer);
                } 
            } while (stepThroughSpacerPath(walk_elem, &previous_spacer));
            
            // if we get to this point then it means that we reached a cross node or the end of a path
            // the first node in the walking elem would be in the current contig
            
            current_contig_spacers.push_back(walk_elem->first());
            
            if ((walk_elem->second())->getSpacerRank() == 1) 
            {
                // end of path
                current_contig_spacers.push_back(walk_elem->second());
            } 
            else 
            {
                // push the cross node onto a list
                cross_nodes.push_back(walk_elem->second());
            }
            
            // assign the nodes the same contig id as the cap -- but not the cross node
            setContigIDForSpacers(&current_contig_spacers);
        } 
        cap_node_iter++;
    }
    NM_NextContigID++;
    walkFromCross(&cross_nodes);

    delete walk_elem;
    
    logInfo("Made: " << NM_NextContigID << " spacer contig(s)", 1);
    return 0;
}

bool NodeManager::walkFromCross(SpacerInstanceList * crossNodes)
{
    WalkingManager * walk_elem = new WalkingManager();
    SpacerInstanceList_Iterator cross_node_iter = crossNodes->begin();
    while (cross_node_iter != crossNodes->end()) 
    {    
        (*cross_node_iter)->setContigID(NM_NextContigID++);
        // walk along those edges as long as possible making sure that the first edge is not a 
        // cross node and that the nodes don't already have a contig id
        SpacerEdgeVector_Iterator edge_iter = (*cross_node_iter)->begin();
        while (edge_iter != (*cross_node_iter)->end()) 
        {
            // go through all the edges of the cross node that do not have a contig id
            if (((*edge_iter)->edge->isAttached()) && (0 == (*edge_iter)->edge->getContigID())) 
            {
                if( getSpacerEdgeFromCross(walk_elem, (*edge_iter)->edge))
                {

                    SpacerInstanceVector current_contig_nodes;
                    SpacerInstance * previous_node = NULL;
                    // edge is a path node so walk
                    do {
                        if (NULL != previous_node) 
                        {
                            current_contig_nodes.push_back(previous_node);
                        }
                    } while (stepThroughSpacerPath(walk_elem, &previous_node));
                    
                    if ((walk_elem->second())->getSpacerRank() == 1 && (walk_elem->second())->isAttached()) 
                    {
                        // end of path
                        current_contig_nodes.push_back(walk_elem->second());
                    } 
                    else if ((walk_elem->second())->getContigID() == 0 && (walk_elem->second())->isAttached())
                    {
                        // add the first node into the list of current contig
                        current_contig_nodes.push_back(walk_elem->first());
                        
                        // push the cross node onto a vector
                        crossNodes->push_back(walk_elem->second());
                    }
                    //current_contig_nodes.push_back(walk_elem->first());
                    setContigIDForSpacers(&current_contig_nodes);
                    NM_NextContigID++;
                }
                else 
                {                    
                    // means that the edge is a cross node so push it back on to the list
                    crossNodes->push_back((*edge_iter)->edge);
                }
            }
            ++edge_iter;
        }
        cross_node_iter++;
    }
    delete walk_elem;
    return true;
    
}


void NodeManager::setContigIDForSpacers(SpacerInstanceVector * currentContigNodes)
{
    SpacerInstanceVector_Iterator iter = currentContigNodes->begin();
    while (iter != currentContigNodes->end()) 
    {
        (*iter)->setContigID(NM_NextContigID);
        iter++;
    }
}

// Printing / IO

void NodeManager::dumpReads(std::string readsFileName, bool showDetached, bool split)
{
    //-----
    // dump reads to this file
    //
    std::set<std::string> reads_set; 
    std::ofstream reads_file;
    reads_file.open(readsFileName.c_str());
    if (reads_file.good()) 
    {
        SpacerListIterator spacer_iter = NM_Spacers.begin();
        while(spacer_iter != NM_Spacers.end())
        {
            SpacerInstance * SI = spacer_iter->second;
            
            // get the crisprnodes for these guys
            CrisprNode * Cleader = SI->getLeader();
            CrisprNode * Clast = SI->getLast();
            
            if(showDetached || ((SI->getLeader())->isAttached() && (SI->getLast())->isAttached()))
            {
                std::vector<std::string> headers = Cleader->getReadHeaders(&NM_StringCheck);
                std::vector<std::string>::iterator h_iter = headers.begin();
                while(h_iter != headers.end())
                {
                    reads_set.insert(*h_iter);
                    h_iter++;
                }
                
                headers = Clast->getReadHeaders(&NM_StringCheck);
                h_iter = headers.begin();
                while(h_iter != headers.end())
                {
                    reads_set.insert(*h_iter);
                    h_iter++;
                }
            }
            spacer_iter++;
        }
        
        // now we can print all the reads to file
        ReadListIterator read_iter = NM_ReadList.begin();
        while (read_iter != NM_ReadList.end()) 
        {
            std::string header = (*read_iter)->getHeader();
            if(reads_set.find(header) != reads_set.end())
            {
                reads_file <<*(*read_iter)<<std::endl;
            }
            read_iter++;
        }
        reads_file.close();
    }
}


// Spacer dictionaries
void NodeManager::addSpacersToDOM(crispr::xml::writer * xmlDoc, 
                                  xercesc::DOMElement * parentNode, 
                                  bool showDetached, 
                                  std::set<StringToken>& allSources)
{
    SpacerListIterator spacer_iter = NM_Spacers.begin();
    while(spacer_iter != NM_Spacers.end())
    {
        SpacerInstance * SI = spacer_iter->second;
        if((showDetached || ((SI->getLeader())->isAttached() && (SI->getLast())->isAttached())) && !(SI->isFlanker()))
        {
            std::set<StringToken> nr_tokens;
            getHeadersForSpacers(SI, nr_tokens);
            
            // generate the spacer tag
            std::string spacer = NM_StringCheck.getString(SI->getID());
            std::string spid = "SP" + to_string(SI->getID());
            std::string cov = to_string(SI->getCount());
            xercesc::DOMElement * spacer_node = xmlDoc->addSpacer(spacer, spid, parentNode, cov);
            appendSourcesForSpacer(spacer_node, 
                                   nr_tokens, 
                                   xmlDoc);
            allSources.insert(nr_tokens.begin(), nr_tokens.end());

        }
        spacer_iter++;
    }
}

void NodeManager::addFlankersToDOM(crispr::xml::writer * xmlDoc, 
                                   xercesc::DOMElement * parentNode, 
                                   bool showDetached, 
                                   std::set<StringToken>& allSources)
{
    SpacerInstanceVector_Iterator iter;
    for (iter = NM_FlankerNodes.begin(); iter != NM_FlankerNodes.end(); iter++) {
        SpacerInstance * SI = *iter;
        if(showDetached || ((SI->getLeader())->isAttached() && (SI->getLast())->isAttached()))
        {
            std::set<StringToken> nr_tokens;
            getHeadersForSpacers(SI, nr_tokens);
            
            std::string spacer = NM_StringCheck.getString(SI->getID());
            std::string flid = "FL" + to_string(SI->getID());
            xercesc::DOMElement * spacer_node = xmlDoc->addFlanker(spacer, flid, parentNode);
            // add in all the source tags for this spacer
            appendSourcesForSpacer(spacer_node, 
                                   nr_tokens, 
                                   xmlDoc);
            allSources.insert(nr_tokens.begin(), nr_tokens.end());

        }
    }
}

void NodeManager::printAssemblyToDOM(crispr::xml::writer * xmlDoc, xercesc::DOMElement * parentNode, bool showDetached)
{
    
    int current_contig_num = 0;
    while (current_contig_num < NM_NextContigID) 
    {
        current_contig_num++;
        std::string cid = "C" + to_string(current_contig_num);
        xercesc::DOMElement * contig_elem = xmlDoc->addContig(cid, parentNode);

        SpacerListIterator spacer_iter = NM_Spacers.begin();
        while(spacer_iter != NM_Spacers.end())
        {
            SpacerInstance * SI = spacer_iter->second;
            if (SI->getContigID() == current_contig_num) 
            {
                if( showDetached || SI->isAttached())
                {

                    std::string spacer = NM_StringCheck.getString(SI->getID());
                    std::string id = (SI->isFlanker()) ? "FL" + to_string(SI->getID()) : "SP" + to_string(SI->getID());
                    
                    xercesc::DOMElement * cspacer = xmlDoc->addSpacerToContig(id, contig_elem);

                    bool ff = false;
                    bool bf = false;
                    bool fs = false;
                    bool bs = false;
                    
                    xercesc::DOMElement * fspacers = NULL;
                    xercesc::DOMElement * bspacers = NULL;
                    xercesc::DOMElement * fflankers = NULL;
                    xercesc::DOMElement * bflankers = NULL;
                    SpacerEdgeVector_Iterator sp_iter = SI->begin();
                    while (sp_iter != SI->end()) 
                    {
                        if ((*sp_iter)->edge->isAttached()) 
                        {

                            std::string edge_id = (SI->isFlanker()) ? "FL" + to_string((*sp_iter)->edge->getID()) : "SP" + to_string((*sp_iter)->edge->getID());
                            std::string drid = "DR1";
                            std::string drconf = "0";
                            std::string directjoin = "0";
                            switch ((*sp_iter)->d) 
                            {
                                case FORWARD:
                                {
                                    if ((*sp_iter)->edge->isFlanker()) {
                                        if (ff) 
                                        {
                                            // we've already created <fflankers>
                                            // add spacer
                                            xmlDoc->addFlanker("ff", edge_id, drconf, directjoin, fflankers);//("fs", edge_spid, drid, drconf, fspacers);
                                        } 
                                        else 
                                        {
                                            // create <fflankers>
                                            fflankers = xmlDoc->createFlankers("fflankers");
                                            xmlDoc->addFlanker("ff", edge_id, drconf, directjoin, fflankers);//("fs", edge_spid, drid, drconf, fspacers);
                                            ff = true;
                                        }
                                    } else {
                                        if (fs) 
                                        {
                                            // we've already created <fspacers>
                                            // add spacer
                                            xmlDoc->addSpacer("fs", edge_id, drid, drconf, fspacers);
                                        } 
                                        else 
                                        {
                                            // create <fspacers>
                                            fspacers = xmlDoc->createSpacers("fspacers");
                                            xmlDoc->addSpacer("fs", edge_id, drid, drconf, fspacers);
                                            fs = true;
                                        }
                                    }

                                    break;
                                }
                                case REVERSE:
                                {
                                    if ((*sp_iter)->edge->isFlanker()) {
                                        if (bf) 
                                        {
                                            // we've already created <bflankers>
                                            // add spacer
                                            xmlDoc->addFlanker("bf", edge_id, drconf, directjoin, bflankers);//("bs", edge_id, drid, drconf, bspacers);
                                        } 
                                        else 
                                        {
                                            // create <bflankers>
                                            bflankers = xmlDoc->createFlankers("bflankers");
                                            xmlDoc->addFlanker("bf", edge_id, drconf, directjoin, bflankers);//("bs", edge_id, drid, drconf, bspacers);
                                            bf = true;
                                        }
                                    } else {
                                        if (bs) 
                                        {
                                            // we've already created <fspacers>
                                            // add spacer
                                            xmlDoc->addSpacer("bs", edge_id, drid, drconf, bspacers);
                                        } 
                                        else 
                                        {
                                            // create <bspacers>
                                            bspacers = xmlDoc->createSpacers("bspacers");
                                            xmlDoc->addSpacer("bs", edge_id, drid, drconf, bspacers);
                                            bs = true;
                                        }
                                    }

                                    break;
                                }
                                default:
                                {
                                    break;
                                }
                            }
                        }
                        ++sp_iter;
                    }
                    if (bspacers != NULL) 
                    {
                        cspacer->appendChild(bspacers);

                    }
                    if (fspacers != NULL) 
                    {
                        cspacer->appendChild(fspacers);

                    }
                    if (bflankers != NULL) 
                    {
                        cspacer->appendChild(bflankers);
                        
                    }
                    if (fflankers != NULL) 
                    {
                        cspacer->appendChild(fflankers);
                        
                    }
                }
            }
            spacer_iter++;
        }
    }

}

void NodeManager::getHeadersForSpacers(SpacerInstance * SI, std::set<StringToken>& nrTokens)
{
    // go through all the string tokens for both the leader and last nodes
    // for all the CrisprNodes in the Spacers
    CrisprNode * first_node = SI->getLeader();
    CrisprNode * second_node = SI->getLast();
    std::vector<StringToken>::iterator header_iter;
    for (header_iter = first_node->beginHeaders(); header_iter != first_node->endHeaders(); header_iter++) {
        nrTokens.insert(*header_iter);
    }
    for (header_iter = second_node->beginHeaders(); header_iter != second_node->endHeaders(); header_iter++) {
        nrTokens.insert(*header_iter);
    }
}

void NodeManager::appendSourcesForSpacer(xercesc::DOMElement * spacerNode, 
                            std::set<StringToken>& nrTokens,
                            crispr::xml::writer * xmlDoc)
{
    // add in all the source tags for this spacer
    std::set<StringToken>::iterator nr_iter;
    for (nr_iter = nrTokens.begin(); nr_iter != nrTokens.end(); nr_iter++) {
        std::string s = NM_StringCheck.getString(*nr_iter);
        std::string sid = "SO";
        sid += to_string(*nr_iter);
        // add the source to both the current spacer
        // and to the total sources list
        //xmlDoc->addSource(s, sid, allSources);
        xmlDoc->addSpacerSource(sid, spacerNode);
    }
}

void NodeManager::generateAllsourceTags(crispr::xml::writer * xmlDoc, 
                           std::set<StringToken>& allSourcesForNM,
                           xercesc::DOMElement * parentNode
                           )
{
    // add in all the source tags for this spacer
    std::set<StringToken>::iterator nr_iter;
    for (nr_iter = allSourcesForNM.begin(); nr_iter != allSourcesForNM.end(); nr_iter++) {
        std::string s = NM_StringCheck.getString(*nr_iter);
        std::string sid = "SO";
        sid += to_string(*nr_iter);
        xmlDoc->addSource(s, sid, parentNode);
    }
}
// Making purdy colours
void NodeManager::setDebugColourLimits(void)
{
    //-----
    // Make the colurs needed for printing the graphviz stuff
    //
    double max_coverage = 0;
    double min_coverage = 10000000;
    NodeListIterator nl_iter = nodeBegin();
    while (nl_iter != nodeEnd()) 
    {
        int coverage = (nl_iter->second)->getCoverage();
        if (coverage > max_coverage) 
        {
            max_coverage = coverage;
        }
        else if(coverage < min_coverage)
        {
            min_coverage = coverage;
        }
        nl_iter++;
    }
    
    NM_DebugRainbow.setType(NM_Opts->graphColourType);
    
    if (NM_Opts->coverageBins != -1) 
    {
        NM_DebugRainbow.setLimits(min_coverage, max_coverage, NM_Opts->coverageBins);
    } 
    else 
    {
        NM_DebugRainbow.setLimits(min_coverage,max_coverage);
    }
}

void NodeManager::setSpacerColourLimits(void)
{
    //-----
    // Make the colurs needed for printing the graphviz stuff
    //
    double max_coverage = 0;
    double min_coverage = 10000000;
    
    SpacerListIterator sp_iter = NM_Spacers.begin();
    while (sp_iter != NM_Spacers.end()) 
    {
        int coverage = (sp_iter->second)->getCount();
        if (coverage > max_coverage) 
        {
            max_coverage = coverage;
        }
        else if(coverage < min_coverage)
        {
            min_coverage = coverage;
        }
        sp_iter++;
    }
    
    NM_SpacerRainbow.setType(NM_Opts->graphColourType);
    if (NM_Opts->coverageBins != -1) 
    {
        NM_SpacerRainbow.setLimits(min_coverage, max_coverage, NM_Opts->coverageBins);
    } 
    else 
    {
        NM_SpacerRainbow.setLimits(min_coverage,max_coverage);
    }
}

void NodeManager::printDebugGraph(std::ostream &dataOut, std::string title, bool showDetached, bool printBackEdges, bool longDesc)
{
    //-----
    // Print a graphviz style graph of the DRs and spacers
    //
    setDebugColourLimits();
    
    gvGraphHeader(dataOut, title);
    NodeListIterator nl_iter = nodeBegin();
    // first loop to print out the nodes
    while (nl_iter != nodeEnd()) 
    {
        // check whether we should print
        if((nl_iter->second)->isAttached() | showDetached)
        {
            printDebugNodeAttributes(dataOut, nl_iter->second ,NM_DebugRainbow.getColour((nl_iter->second)->getCoverage()), longDesc);
        }
        nl_iter++;
    }
    
    // and go through again to print the edges
    nl_iter = nodeBegin();
    while (nl_iter != nodeEnd()) 
    {
        // check whether we should print
        if((nl_iter->second)->isAttached() | showDetached)
        {
            std::stringstream ss;
            if(longDesc)
                ss << (nl_iter->second)->getID() << "_" << NM_StringCheck.getString((nl_iter->second)->getID());
            else
                ss << (nl_iter->second)->getID();
            (nl_iter->second)->printEdges(dataOut, &NM_StringCheck, ss.str(), showDetached, printBackEdges, longDesc);
        }
        nl_iter++;
    }
    gvGraphFooter(dataOut)
}

void NodeManager::printDebugNodeAttributes(std::ostream& dataOut, CrisprNode * currCrisprNode, std::string colourCode, bool longDesc)
{
    //-----
    // print the node declaration
    //
    std::stringstream ss;
    if(longDesc)
        ss << currCrisprNode->getID() << "_" << NM_StringCheck.getString(currCrisprNode->getID());
    else
        ss << currCrisprNode->getID();
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

bool NodeManager::printSpacerGraph(std::string& outFileName, std::string title, bool longDesc, bool showSingles)
{
    //-----
    // Print a graphviz style graph of the DRs and spacers
    //
    std::stringstream tmp_out;
    setSpacerColourLimits();
    gvGraphHeader(tmp_out, title);
    bool at_least_one_spacer=false;        
    SpacerListIterator spi_iter = NM_Spacers.begin();
    while (spi_iter != NM_Spacers.end()) 
    {
        if ((spi_iter->second)->isAttached() && (showSingles || (0 != (spi_iter->second)->getSpacerRank()))) 
        {
            at_least_one_spacer = true;
            // print the graphviz nodes
            std::string label = getSpacerGraphLabel(spi_iter->second, longDesc);
            
            // print the node attribute
            if ((spi_iter->second)->isFlanker()) {
                gvFlanker(tmp_out, label, NM_SpacerRainbow.getColour((spi_iter->second)->getCount()));
            } else {
                gvSpacer(tmp_out,label,NM_SpacerRainbow.getColour((spi_iter->second)->getCount()));
                
            }
        }
        spi_iter++;
    }
    if (!at_least_one_spacer) 
    {
        return false;
    }  
    
    std::ofstream data_out;
    data_out.open(outFileName.c_str());
    if (data_out.good()) 
    {
        data_out<<tmp_out.str();
        spi_iter = NM_Spacers.begin();
        while (spi_iter != NM_Spacers.end()) 
        {
            if ((spi_iter->second)->isAttached() && (showSingles || (0 != (spi_iter->second)->getSpacerRank()))) 
            {
                // print the graphviz nodes
                std::string label = getSpacerGraphLabel(spi_iter->second, longDesc);
                // print the node attribute
                // now print the edges
                SpacerEdgeVector_Iterator edge_iter = (spi_iter->second)->begin();
                while (edge_iter != (spi_iter->second)->end()) 
                {
                    if (((*edge_iter)->edge)->isAttached() && (*edge_iter)->d == FORWARD && (showSingles || (0 != ((*edge_iter)->edge)->getSpacerRank()))) 
                    {
                        
                        // get the label for our edge
                        // print the graphviz nodes
                        gvSpEdge(data_out, label, getSpacerGraphLabel((*edge_iter)->edge, longDesc));
                    }
                    edge_iter++;
                }
            }   
            spi_iter++;
        }        
        gvGraphFooter(data_out);
        data_out.close();
        return true;
    } 
    else 
    {
        logError("Cannot open output file "<<outFileName);
        return false;
    }
}
    
std::string NodeManager::getSpacerGraphLabel(SpacerInstance * spacer, bool longDesc)
{
    //-----
    // Get the label for a spacer when printing the spacer graph
    //
    std::stringstream se;
    if(longDesc)
    {
        if (spacer->isFlanker()) {
            se<< CRASS_DEF_GV_FL_PREFIX;
        } else {
            se << CRASS_DEF_GV_SPA_PREFIX;
        }
        se << spacer->getID() << "_" << NM_StringCheck.getString(spacer->getID()) << "_" << spacer->getCount();
    }
    else
    {
        if (spacer->isFlanker()) {
            se<< CRASS_DEF_GV_FL_PREFIX;
        } else {
            se << CRASS_DEF_GV_SPA_PREFIX;
        }
        se << spacer->getID() << "_" << spacer->getCount();
    }
    se << "_C" << spacer->getContigID();
    return se.str();
}

void NodeManager::printAllSpacers(void)
{
    SpacerListIterator sp_iter = NM_Spacers.begin();
    
    // first we build up the contents of the walking queue
    while(sp_iter != NM_Spacers.end())
    {
        (sp_iter->second)->printContents(); 
        sp_iter++;
    }
}

void NodeManager::printSpacerKey(std::ostream &dataOut, int numSteps, std::string groupNumber)
{
    //-----
    // Print a graphviz style graph of the DRs and spacers
    //
    static int cluster_number = 0;
    gvKeyGroupHeader(dataOut, cluster_number, groupNumber);
    double ul = NM_SpacerRainbow.getUpperLimit();
    double ll = NM_SpacerRainbow.getLowerLimit();
    double step_size = (ul - ll) / (numSteps - 1);
    if(step_size < 1) { step_size = 1; }
    for(double i = ll; i <= ul; i+= step_size)
    {
        int this_step = int(i); 
        std::stringstream ss;
        ss << this_step;
        gvKeyEntry(dataOut, ss.str(), NM_SpacerRainbow.getColour(this_step));
    }
    gvKeyFooter(dataOut);
    cluster_number++;
}

// Flankers

void NodeManager::generateFlankers(bool showDetached)
{
    // generate statistics about length of spacers
    getSpacerCountAndStats();
    
    
    // do some maths
    double stdev = NM_SpacerLenStat.standardDeviation();
    int mean = (int)NM_SpacerLenStat.mean();
    int lower_bound = mean - (stdev*1.5);
    int upper_bound = mean + (stdev*1.5);
    
    // if there is no variation then there is no flankers
    if (stdev > 1 ) 
    {
        logInfo("Ave SP Length: "<<mean<<" Deviation: "<<stdev<<" UB: "<<upper_bound<<" LB: "<<lower_bound, 3);
        // call a spacer a 'flanker' if it's length is more than 1 standard deviation from the mean length
        // and it is a cap node

        SpacerListIterator spacer_iter = NM_Spacers.begin();
        //spacer_iter = NM_Spacers.begin();
        while(spacer_iter != NM_Spacers.end())
        {
            SpacerInstance * SI = spacer_iter->second;
 
            if(showDetached || ((SI->getLeader())->isAttached() && (SI->getLast())->isAttached()))
            {
                /*if(SI->isCap()) {*/
                    int spacer_length = (int)(NM_StringCheck.getString(SI->getID())).length();
                    if (spacer_length > upper_bound || spacer_length < lower_bound) {
                        SI->setFlanker(true);
                        NM_FlankerNodes.push_back(SI);
                   }
                /* }*/
            }
            spacer_iter++;
        }
    }
    else 
    {
        logInfo("not enough length variation to detect flankers", 3);
    }
    
    // do this here so that any subsuquent calls don't include the flanking sequences
    clearStats();
}
