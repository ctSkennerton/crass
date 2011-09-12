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

NodeManager::NodeManager(std::string drSeq, StringCheck * strCheck)
{
    //-----
    // constructor
    //
    NM_DirectRepeatSequence = drSeq;
    NM_StringCheck = strCheck;
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

    if(RH->getFirstSpacer(&working_str))
    {
        do {
           addCrisprNodes(&prev_node, working_str);
        } while(RH->getNextSpacer(&working_str));
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
void NodeManager::addCrisprNodes(CrisprNode ** prevNode, std::string& workingString)
{
    //-----
    // Given a spacer string, cut kmers from each end and make crispr nodes
    //
    // now cut kmers on either side and add the pair into the node manager 
    std::string first_kmer = workingString.substr(0,CRASS_DEF_NODE_KMER_SIZE);
    std::string second_kmer = workingString.substr(workingString.length() - CRASS_DEF_NODE_KMER_SIZE, CRASS_DEF_NODE_KMER_SIZE );
    CrisprNode * first_kmer_node;
    CrisprNode * second_kmer_node;
    
    // we need to know if we've seen both of the guys before
    bool seen_first = false;
    bool seen_second = false;

    // check to see if these kmers are already stored
    StringToken st1 = NM_StringCheck->getToken(first_kmer);
    StringToken st2 = NM_StringCheck->getToken(second_kmer);

    // if they have been added previously then token != 0
    if(0 == st1)
    {
        // first time we've seen this guy. Make some new objects
        st1 = NM_StringCheck->addString(first_kmer);
        first_kmer_node = new CrisprNode(st1);
        
        // add them to the pile
        NM_Nodes[st1] = first_kmer_node;
    }
    else
    {
        // we already have a node for this guy
        first_kmer_node = NM_Nodes[st1];
        (NM_Nodes[st1])->incrementCount();
        
        seen_first = true;
    }
    if(0 == st2)
    {
        st2 = NM_StringCheck->addString(second_kmer);
        second_kmer_node = new CrisprNode(st2);
        second_kmer_node->setForward(false);
        NM_Nodes[st2] = second_kmer_node;
    }
    else
    {
        second_kmer_node = NM_Nodes[st2];
        (NM_Nodes[st2])->incrementCount();

        seen_second = true;
    }
    
    // the first kmers pair is the previous node which lay before it therefore bool is true
    // make sure prevNode is not NULL
    if (*prevNode != NULL && !seen_first) 
    {
        (*prevNode)->addEdge(first_kmer_node, CN_EDGE_JUMPING_F);
        first_kmer_node->addEdge(*prevNode, CN_EDGE_JUMPING_B);
    }

    if(!(seen_first & seen_second))
    {
        first_kmer_node->addEdge(second_kmer_node, CN_EDGE_FORWARD);
        second_kmer_node->addEdge(first_kmer_node, CN_EDGE_BACKWARD);
    }
    
    // now it's time to add the spacer
    SpacerInstance * curr_spacer;

    // check to see if we already have it here
    SpacerKey this_sp_key = makeSpacerKey(st1, st2);

    if(NM_Spacers.find(this_sp_key) == NM_Spacers.end())
    {
        // new instance
        StringToken sp_str_token = NM_StringCheck->addString(workingString);
        curr_spacer = new SpacerInstance(sp_str_token, first_kmer_node, second_kmer_node);
        NM_Spacers[this_sp_key] = curr_spacer;
    }
    else
    {
        // increment the number of times we've seen this guy
        (NM_Spacers[this_sp_key])->incrementCount();
    }

    *prevNode = second_kmer_node;
}

// Walking


// Cleaning
void NodeManager::cleanGraph(void)
{
    //-----
    // Clean all the bits off the graph mofo!
    //
    NodeListIterator bob = nodeBegin();
    bob++;
    bob++;
    bob++;
    bob++;
    bob++;
    bob++;
    (bob->second)->detachNode();
}

// Making purdy colours
void NodeManager::setUpperAndLowerCoverage(void)
{
    // loop through all of the nodes and determine the upper and lower dounds for our graph
    NodeListIterator nl_iter = nodeBegin();
    while (nl_iter != nodeEnd()) 
    {
        int coverage = (nl_iter->second)->getCoverage();
        if (coverage > NM_MaxCoverage) 
        {
            NM_MaxCoverage = coverage;
        }
        else if(coverage < NM_MinCoverage)
        {
            NM_MinCoverage = coverage;
        }
        nl_iter++;
    }
}


void NodeManager::setColourLimits(void)
{
    //-----
    // Make the colurs needed for printing the graphviz stuff
    //
    setUpperAndLowerCoverage();
    NM_Rainbow.setLimits(NM_MinCoverage,NM_MaxCoverage);
}

// Printing / IO
void NodeManager::printGraph(std::ostream &dataOut, std::string title, bool showDetached, bool printBackEdges)
{
    //-----
    // Print a graphviz style graph of the DRs and spacers
    //
    gvGraphHeader(dataOut, title);
    NodeListIterator nl_iter = nodeBegin();
    while (nl_iter != nodeEnd()) 
    {
        // check whether we should print
        if((nl_iter->second)->isAttached() || showDetached)
        {
            (nl_iter->second)->printEdges(dataOut, showDetached, printBackEdges, NM_Rainbow.getColour((nl_iter->second)->getCoverage()));
        }
        nl_iter++;
    }
    gvGraphFooter(dataOut)
}
  



















