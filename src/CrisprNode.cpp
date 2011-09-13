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


//
// Edge level functions
//
bool CrisprNode::addEdge(CrisprNode * parterNode, EDGE_TYPE type)
{
    //-----
    // Add a new edge return success if the partner has been added
    // 
    
    //logInfo("Adding "<<sayEdgeTypeLikeAHuman(type)<<" edge from "<<getID()<<" to "<<parterNode->getID(), 1);

    // get the right edge list
    edgeList * add_list = getEdges(type);
    
    // now see we haven't added it before
    if(add_list->find(parterNode) == add_list->end())
    {
        // new guy
        (*add_list)[parterNode] = true;
        return true;
    }
    else
    {
        // already got this guy
        logError("Adding edge ("<<parterNode->getID()<<") twice to CN ("<<getID()<<")");
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
            return &CN_ForwardEdges;
        case CN_EDGE_BACKWARD:
            return &CN_BackwardEdges;
        case CN_EDGE_JUMPING_F:
            return &CN_JumpingForwardEdges;
        case CN_EDGE_JUMPING_B:
            return &CN_JumpingBackwardEdges;
        default:
            logError("Edge type not known! Returning NULL...");
            return NULL;
    }
}

//
// Node level functions
//
void CrisprNode::setAttach(bool attachState)
{
    //-----
    // detach or re-attach this node
    //
    
    // find and attached nodes and set the edges to attachState
    edgeListIterator eli = CN_ForwardEdges.begin();
    while(eli != CN_ForwardEdges.end())
    {
        // go through each edge, check if it's not the right state
        if((eli->second ^ attachState) && (eli->first)->isAttached())
        {
            // this edge is not the right state and the corresponding node is actually attached
            edgeList * other_eli = (eli->first)->getEdges(CN_EDGE_BACKWARD);
            (*other_eli)[this] = attachState;
            eli->second = attachState;
        }
        eli++;
    }
    eli = CN_BackwardEdges.begin();
    while(eli != CN_BackwardEdges.end())
    {
        if((eli->second ^ attachState) && (eli->first)->isAttached())
        {
            edgeList * other_eli = (eli->first)->getEdges(CN_EDGE_FORWARD);
            (*other_eli)[this] = attachState;
            eli->second = attachState;
        }
        eli++;
    }
    eli = CN_JumpingForwardEdges.begin();
    while(eli != CN_JumpingForwardEdges.end())
    {
        if((eli->second ^ attachState) && (eli->first)->isAttached())
        {
            edgeList * other_eli = (eli->first)->getEdges(CN_EDGE_JUMPING_B);
            (*other_eli)[this] = attachState;
            eli->second = attachState;
        }
        eli++;
    }
    eli = CN_JumpingBackwardEdges.begin();
    while(eli != CN_JumpingBackwardEdges.end())
    {
        std::cout<<(eli->first)->isAttached()<<" : "<<eli->second<< " : "<< attachState<<std::endl;
        if((eli->second ^ attachState) && (eli->first)->isAttached())
        {
            edgeList * other_eli = (eli->first)->getEdges(CN_EDGE_JUMPING_F);
            (*other_eli)[this] = attachState;
            eli->second = attachState;
        }
        eli++;
    }
    
    // set our state
    CN_Attached =  attachState;       
}

int CrisprNode::getRank(EDGE_TYPE type)
{
    //-----
    // return the rank of the node
    //
    switch(type)
    {
        case CN_EDGE_FORWARD:
            return CN_InnerRank_F;
        case CN_EDGE_BACKWARD:
            return CN_InnerRank_B;
        case CN_EDGE_JUMPING_F:
            return CN_JumpingRank_F;
        case CN_EDGE_JUMPING_B:
            return CN_JumpingRank_B;
        default:
            logError("Edge type not know! Returning -1...")
            return -1;
    }
}

//
// File IO / printing
//


void CrisprNode::printEdges(std::ostream &dataOut, bool showDetached, bool printBackEdges)
{
    //-----
    // print the edges so that the first member of the pair is first
    //
        
    // now print the edges
    edgeListIterator eli = CN_ForwardEdges.begin();
    while(eli != CN_ForwardEdges.end())
    {
        // check if the edge is active
        if((eli->second) || showDetached)
            gvEdge(dataOut,CN_id,(eli->first)->getID());
        eli++;
    }
    eli = CN_JumpingForwardEdges.begin();
    while(eli != CN_JumpingForwardEdges.end())
    {
         if((eli->second) || showDetached)
            gvJumpingEdge(dataOut,CN_id,(eli->first)->getID());
        eli++;
    }
    if(printBackEdges)
    {
        eli = CN_BackwardEdges.begin();
        while(eli != CN_BackwardEdges.end())
        {
            if((eli->second) || showDetached)
                gvEdge(dataOut,CN_id,(eli->first)->getID());
            eli++;
        }
        eli = CN_JumpingBackwardEdges.begin();
        while(eli != CN_JumpingBackwardEdges.end())
        {
            if((eli->second) || showDetached)
                gvJumpingEdge(dataOut,CN_id,(eli->first)->getID());
            eli++;
        }
    }
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
            logError("Edge type not known! Returning '---'...");
            return "---";
    }
}


