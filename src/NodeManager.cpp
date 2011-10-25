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
    
    CrisprNode * first_kmer_node;
    CrisprNode * second_kmer_node;
    SpacerKey this_sp_key;
    
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
    if (*prevNode != NULL) 
    {
		this_sp_key = makeSpacerKey(st1, (*prevNode)->getID());
		if(mSpacers.find(this_sp_key) == mSpacers.end())
		{
	        (*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
	        first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
		}
    }

    // now it's time to add the spacer
    SpacerInstance * curr_spacer;

    // check to see if we already have it here
    this_sp_key = makeSpacerKey(st1, st2);

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
    StringToken st2 = mStringCheck.getToken(second_kmer);
    
    // if they have been added previously then token != 0
    if(0 == st2)
    {
        // first time we've seen this guy. Make some new objects
        st2 = mStringCheck.addString(second_kmer);
        second_kmer_node = new CrisprNode(st2);
        second_kmer_node->setForward(false);

        // add them to the pile
        mNodes[st2] = second_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        second_kmer_node = mNodes[st2];
        (mNodes[st2])->incrementCount();
    }
    
    // add in the read headers for the this CrisprNode
    second_kmer_node->addReadHeader(headerSt);
    
    // add this guy in as the previous node for the next iteration
    *prevNode = second_kmer_node;
    
    // there is no one yet to make an edge
}

void NodeManager::addFirstCrisprNode(CrisprNode ** prevNode, std::string& workingString, StringToken headerSt)
{
    if (workingString.length() < CRASS_DEF_NODE_KMER_SIZE)
        return;

    std::string first_kmer = workingString.substr(0, CRASS_DEF_NODE_KMER_SIZE);
    CrisprNode * first_kmer_node;
    
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
	//-----
	// find cap nodes!
	//
	capNodes->clear();

    NodeListIterator all_node_iter = mNodes.begin();
    while (all_node_iter != mNodes.end()) 
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

	NodeListIterator all_node_iter = mNodes.begin();
	while (all_node_iter != mNodes.end()) 
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
	
    NodeListIterator all_node_iter = mNodes.begin();
    while (all_node_iter != mNodes.end()) 
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

// Cleaning
int NodeManager::cleanGraph(void)
{
    //-----
    // Clean all the bits off the graph mofo!
    //
	std::multimap<CrisprNode *, CrisprNode *> fork_choice_map;
	NodeVector nv_cap, nv_other, detach_list;
	NodeVectorIterator nv_iter;
	
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
	
	// then do bubbles
	nv_iter = nv_other.begin();
	while(nv_iter != nv_other.end())
	{
		// get the total rank first.  A linear node should have 4 edges
		if((*nv_iter)->getTotalRank() != 4)
		{
			// get the rank for the the inner and jumping edges.
			if((*nv_iter)->getInnerRank() != 2)
			{
				// there are multiple inner edges for this guy
                clearBubbles(*nv_iter, CN_EDGE_FORWARD);

			}
			else if((*nv_iter)->getJumpingRank() != 2)
			{
				// there are multiple jumping edges for this guy
                clearBubbles(*nv_iter, CN_EDGE_JUMPING_F);
			}
			else
			{
				// TODO: something has gone wrong
			}
		}
		
		nv_iter++;
	}
	
	// finally, detach!
	nv_iter = detach_list.begin();
	while(nv_iter != detach_list.end())
	{
		(*nv_iter)->detachNode();
		nv_iter++;
	}

	return 0;
}

void NodeManager::clearBubbles(CrisprNode * rootNode, EDGE_TYPE currentEdgeType)
{
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
            int key = makeKey(rootNode->getID(),(curr_edges_iter->first)->getID());
            
            bubble_map[key] = (curr_edges_iter->first)->getID();
            
            // we also want to go through all the edges of the nodes above (2nd degree separation)
            // and since we used the forward edges to get here we now want the opposite (Jummping_F)
            edgeList * edges_of_curr_edge = (curr_edges_iter->first)->getEdges(currentEdgeType);
            
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
                        bubble_map[new_key] = (edges_of_curr_edge_iter->first)->getID();
                    } 
                    else 
                    {
                        // aha! he is pointing back onto the same guy as someone else.  We have a bubble!
                        //get the CrisprNode of the first guy
                        
                        CrisprNode * first_node = mNodes[bubble_map[new_key]];
                        
                        
#ifdef DEBUG
                        logInfo("Bubble found conecting "<<rootNode->getID()<<" : "<<first_node->getID()<<" : "<<(edges_of_curr_edge_iter->first)->getID(), 7);
#endif
                        //perform a coverage test on the nodes that end up here and kill the one with the least coverage
                        
                        // TODO: this is a pretty dumb test, since the coverage between the two nodes could be very similar
                        // for example one has a coverage of 20 and the other a coverage of 18.  The node with 18 will get
                        // removed but should it have been?  The way to fix this would be to create some infastructure in
                        // NodeManager to calculate the average and stdev of the coverage and then remove a node only if
                        // it is below 1 stdev of the average, else it could be a biological thing that this bubble exists.
                        
                        if (first_node->getCoverage() > (edges_of_curr_edge_iter->first)->getCoverage()) 
                        {
                            // the first guy has greater coverage so detach our current node
                            (edges_of_curr_edge_iter->first)->detachNode();
#ifdef DEBUG
                            logInfo("Detaching "<<(edges_of_curr_edge_iter->first)->getID()<<" as it has lower coverage", 8);
#endif
                            //TODO: and update the bool in the edgelist?
                            // edges_of_curr_edge_iter->second = false;
                        } 
                        else 
                        {
                            // the first guy was lower so kill him
                            first_node->detachNode();
#ifdef DEBUG
                            logInfo("Detaching "<<first_node->getID()<<" as it has lower coverage", 8);
#endif

                            //TODO: need to find him in the edge list in the outer while loop and set to false?
                            // curr_edges->at(first_node) = false;
                        }
                    }
                }
                edges_of_curr_edge_iter++;
            }
        }
        curr_edges_iter++;
    }
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
#ifdef DEBUG
    logInfo("Max Node Coverage: "<<mMaxCoverage<<" Min Node Coverage: "<<mMinCoverage<<std::endl,5);
#endif
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


void NodeManager::printGraph(std::ostream &dataOut, std::string title, bool showDetached, bool printBackEdges, bool longDesc)
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
            printNodeAttributes(dataOut, nl_iter->second ,mRainbow.getColour((nl_iter->second)->getCoverage()), longDesc);
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
        	if(longDesc)
        		ss << (nl_iter->second)->getID() << "_" << mStringCheck.getString((nl_iter->second)->getID());
        	else
        		ss << (nl_iter->second)->getID();
            (nl_iter->second)->printEdges(dataOut, &mStringCheck, ss.str(), showDetached, printBackEdges, longDesc);
        }
        nl_iter++;
    }
    gvGraphFooter(dataOut)
}

void NodeManager::printNodeAttributes(std::ostream& dataOut, CrisprNode * currCrisprNode, std::string colourCode, bool longDesc)
{
	//-----
    // print the node declaration
	//
	std::stringstream ss;
	if(longDesc)
		ss << currCrisprNode->getID() << "_" << mStringCheck.getString(currCrisprNode->getID());
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

int NodeManager::printReadInfo(void)
{
    //-----
    // Print information about the reads from this group
    //
	std::multimap<std::string, StringToken> seen_map;
	NodeVector nv_cap, nv_other;
	
	// get the cap nodes
	findAllNodes(&nv_cap, &nv_other);
	
	// go through once and build the multimap
	NodeVectorIterator nv_iter = nv_cap.begin();
	while(nv_iter != nv_cap.end())
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
	nv_iter = nv_cap.begin();
	while(nv_iter != nv_cap.end())
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

















