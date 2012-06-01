// File: CrisprNode.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Individual nodes for each kmer!
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
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

// local includes
#include "CrisprNode.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "GraphDrawingDefines.h"
#include "Rainbow.h"
#include "StringCheck.h"
#include "libcrispr.h"
#include "ReadHolder.h"
#include <libcrispr/Exception.h>

//
// Edge level functions
//
bool CrisprNode::addEdge(CrisprNode * parterNode, EDGE_TYPE type)
{
    //-----
    // Add a new edge return success if the partner has been added
    // 
    
    // get the right edge list
    edgeList * add_list = getEdges(type);
    
    // now see we haven't added it before
    if(add_list->find(parterNode) == add_list->end())
    {
        // new guy
        (*add_list)[parterNode] = true;
        switch(type)
        {
            case CN_EDGE_FORWARD:
                mInnerRank_F++;
                break;
            case CN_EDGE_BACKWARD:
                mInnerRank_B++;
                break;
            case CN_EDGE_JUMPING_F:
                mJumpingRank_F++;
                break;
            case CN_EDGE_JUMPING_B:
                mJumpingRank_B++;
                break;
            default:
            	throw crispr::runtime_exception(__FILE__,
                                                __LINE__,
                                                __PRETTY_FUNCTION__,
                                                "Unknown edge type");
        }
        return true;
    }
    return false;
}

edgeList * CrisprNode::getEdges(EDGE_TYPE type)
{
    //-----
    // get edges of a particular type
    //
    switch(type)
    {
        case CN_EDGE_FORWARD:
            return &mForwardEdges;
        case CN_EDGE_BACKWARD:
            return &mBackwardEdges;
        case CN_EDGE_JUMPING_F:
            return &mJumpingForwardEdges;
        case CN_EDGE_JUMPING_B:
            return &mJumpingBackwardEdges;
        default:
            throw crispr::runtime_exception(__FILE__,
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            "Unknown edge type");
    }
}

int CrisprNode::calculateReadCoverage(edgeList * currentList, std::map<StringToken, int>& countingMap)
{
    int total_inners = 0;
    edgeListIterator eli;
    for (eli = currentList->begin(); eli != currentList->end(); eli++)
    {
    	// check if he's attached
    	if(! eli->second)
    	{
            continue;
        }
        total_inners++;
#ifdef DEBUG
        logInfo("Edge: "<<(eli->first)->getID(), 10);
#endif
        // get the headers
        std::vector<StringToken> * inner_headers = (eli->first)->getReadHeaders();
        std::vector<StringToken>::iterator inner_rh_iter = inner_headers->begin();
        std::vector<StringToken>::iterator inner_rh_last = inner_headers->end();
        while(inner_rh_iter != inner_rh_last)
        {
#ifdef DEBUG
            logInfo("\tHeader: "<<*inner_rh_iter, 10);
#endif
            if(countingMap.find(*inner_rh_iter) != countingMap.end()) {
                countingMap[*inner_rh_iter]++;
#ifdef DEBUG
                logInfo("\t\tIncrementing: "<<*inner_rh_iter<<" : "<<countingMap[*inner_rh_iter], 10);
#endif
            }
            inner_rh_iter++;
        }
    }
    return total_inners;
}


int CrisprNode::getDiscountedCoverage(void)
{
	//-----
	// Return a (possibly) lower version of the coverage
	// used when removing bubbles
	//
	std::map<StringToken, int> counting_map;
	
	// initialise the counting map to inlcude all the reads we care about
	std::vector<StringToken>::iterator rh_iter = mReadHeaders.begin();
	for (rh_iter = mReadHeaders.begin(); rh_iter != mReadHeaders.end(); ++rh_iter) {
		counting_map[*rh_iter] = 0;
#ifdef DEBUG
        logInfo("Node :"<<mid<<" Header: "<<*rh_iter, 10);
#endif
	}
#ifdef DEBUG
    logInfo("Node: "<<mid<<" Headers size:"<<mReadHeaders.size(), 10);
    logInfo("\tForward: "<<mForwardEdges.size(), 10);
    logInfo("\tBackward: "<<mBackwardEdges.size(), 10);
    logInfo("\tJForward: "<<mJumpingForwardEdges.size(), 10);
    logInfo("\tJBackward: "<<mJumpingBackwardEdges.size(), 10);
#endif
	// now update the counting map with reads found on the innner connecting nodes -> perhaps one of these lists is empty?
	int total_inners = 0;
	// first forward
	total_inners += calculateReadCoverage(&mForwardEdges, counting_map);
#ifdef DEBUG
	logInfo("Node "<<mid<<": After forward inners: "<<total_inners, 10);
#endif
	// then backward
	total_inners += calculateReadCoverage(&mBackwardEdges, counting_map);	
#ifdef DEBUG
	logInfo("Node "<<mid<<": After Backward inners: "<<total_inners, 10);
#endif
    int ret_val = 0;
    std::map<StringToken, int>::iterator cm_iter = counting_map.begin();
    std::map<StringToken, int>::iterator cm_last = counting_map.end();
    while(cm_iter != cm_last)
    {
    	if(cm_iter->second == total_inners)
    		ret_val++;
    	cm_iter++;
    }
    
    return ret_val;
}


//
// Node level functions
//
void CrisprNode::setEdgeAttachState(edgeList * currentList, bool attachState, EDGE_TYPE currentType)
{
    edgeListIterator eli;
    for (eli = currentList->begin(); eli != currentList->end(); eli++) {

        // go through each edge, check if it's not the right state
        if((eli->second ^ attachState) && (eli->first)->isAttached())
        {
            // this edge is not the right state and the corresponding node is actually attached
            edgeList * other_eli = (eli->first)->getEdges(currentType);
            (*other_eli)[this] = attachState;
            eli->second = attachState;
            (eli->first)->updateRank(attachState, currentType);
            if((eli->first)->getTotalRank() == 0)
            	(eli->first)->setAsDetached();
        }
    }
}

void CrisprNode::setAttach(bool attachState)
{
    //-----
    // detach or re-attach this node
    //
    
    // find and attached nodes and set the edges to attachState
    setEdgeAttachState(&mForwardEdges, attachState, CN_EDGE_FORWARD);

    setEdgeAttachState(&mBackwardEdges, attachState, CN_EDGE_BACKWARD);

    setEdgeAttachState(&mJumpingForwardEdges, attachState, CN_EDGE_JUMPING_F);

    setEdgeAttachState(&mJumpingBackwardEdges, attachState, CN_EDGE_JUMPING_B);
    
    // set our state
    mAttached =  attachState;       
}

int CrisprNode::getRank(EDGE_TYPE type)
{
    //-----
    // return the rank of the node
    //
    switch(type)
    {
        case CN_EDGE_FORWARD:
            return mInnerRank_F;
        case CN_EDGE_BACKWARD:
            return mInnerRank_B;
        case CN_EDGE_JUMPING_F:
            return mJumpingRank_F;
        case CN_EDGE_JUMPING_B:
            return mJumpingRank_B;
        default:
            throw crispr::runtime_exception(__FILE__,
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            "Unknown edge type");
    }
}

void CrisprNode::updateRank(bool attachState, EDGE_TYPE type)
{
    //-----
	// increment or decrement the rank of this type
    //
	int increment = -1;
	if(attachState)
		increment = 1;
	switch(type)
	{
		case CN_EDGE_FORWARD:
			mInnerRank_F += increment;
            break;
		case CN_EDGE_BACKWARD:
			mInnerRank_B += increment;
            break;
		case CN_EDGE_JUMPING_F:
			mJumpingRank_F += increment;
            break;
		case CN_EDGE_JUMPING_B:
			mJumpingRank_B += increment;
            break;
        default:
        	throw crispr::runtime_exception(__FILE__,
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            "Unknown edge type");
	}
}

//
// File IO / printing
//
void CrisprNode::printEdgesForList(edgeList * currentList,
                       std::ostream &dataOut,
                       StringCheck * ST,
                       std::string label, 
                       bool showDetached, 
                       bool longDesc)
{
    edgeListIterator eli; 
    for (eli = currentList->begin(); eli != currentList->end(); eli++) {
        // check if the edge is active
        if((eli->second) || showDetached)
        {
        	std::stringstream ss;
        	if(longDesc)
        		ss << (eli->first)->getID() << "_" << ST->getString((eli->first)->getID());
        	else
        		ss << (eli->first)->getID();
            gvEdge(dataOut,label,ss.str());
        }
    }
}



void CrisprNode::printEdges(std::ostream &dataOut, 
                            StringCheck * ST, 
                            std::string label, 
                            bool showDetached, 
                            bool printBackEdges, 
                            bool longDesc)
{
    //-----
    // print the edges so that the first member of the pair is first
    //
        
    // now print the edges
    printEdgesForList(&mForwardEdges, dataOut, ST, label, showDetached, longDesc);
    printEdgesForList(&mJumpingForwardEdges, dataOut, ST, label, showDetached, longDesc);
    
    if(printBackEdges)
    {
        printEdgesForList(&mBackwardEdges, dataOut, ST, label, showDetached, longDesc);

        printEdgesForList(&mJumpingBackwardEdges, dataOut, ST, label, showDetached, longDesc);

    }
}

std::vector<std::string> CrisprNode::getReadHeaders(StringCheck * ST)
{
	//-----
	// Given a StringCheck return the read headers for this node
	//
	std::vector<std::string> ret_vector;
	if(ST != NULL)
	{
		std::vector<StringToken>::iterator rh_iter = mReadHeaders.begin();
		while(rh_iter != mReadHeaders.end())
		{
			ret_vector.push_back(ST->getString(*rh_iter));
			rh_iter++;
		}
	}
	return ret_vector;
}

std::string CrisprNode::sayEdgeTypeLikeAHuman(EDGE_TYPE type)
{
    //-----
    // get edges of a particular type
    //
    switch(type)
    {
        case CN_EDGE_FORWARD:
            return "Forward";
        case CN_EDGE_BACKWARD:
            return "Backward";
        case CN_EDGE_JUMPING_F:
            return "JumpingForward";
        case CN_EDGE_JUMPING_B:
            return "JumpingBackward";
        default:
            throw crispr::runtime_exception(__FILE__,
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            "Unknown edge type");
    }
}


