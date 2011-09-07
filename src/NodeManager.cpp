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
#include <string>
#include <sstream>
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
    
    // clean up al lthe cripsr nodes
    NodeListIterator node_iter = NM_Nodes.begin();
    while(node_iter != NM_Nodes.end())
    {
        if(NULL != *node_iter)
        {
            delete *node_iter;
            *node_iter = NULL;
        }
        node_iter++;
    }
    NM_Nodes.clear();
    
    SpacerListIterator spacer_iter = NM_Spacers.begin();
    while(spacer_iter != NM_Spacers.end())
    {
        if(NULL != *spacer_iter)
        {
            delete *spacer_iter;
            *spacer_iter = NULL;
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
    this->splitReadHolder(RH);
    NM_ReadList.push_back(RH);
    return true;
}

//----
// private function called from addReadHolder to split the read into spacers and pass it through to others
//
void NodeManager::splitReadHolder(ReadHolder * RH)
{

    std::string working_str;
    CrisprNode * prev_node = NULL;
    
    int size_of_RH_startStops = RH->drListSize();
    if(0 == RH->front())
    {
        if(RH->getFirstSpacer(&working_str))
        {
            addCrisprNodes(&prev_node, working_str);
            size_of_RH_startStops--;
        }
        else
        {
            logError("could't get the first spacer");
        }
        
        while(size_of_RH_startStops > 0)
        {

            if(RH->getNextSpacer(&working_str))
            {
                addCrisprNodes(&prev_node, working_str);
                size_of_RH_startStops -= 2;
            }
            else
            {
                break;
            }
        }
        
        // do the last spacer if required
        // TODO  -1?
        if (RH->back() <= RH->seqLength() - CRASS_DEF_NODE_KMER_SIZE) 
        {
            RH->getNextSpacer(&working_str);
            
            std::string first_kmer = working_str.substr(0, CRASS_DEF_NODE_KMER_SIZE);
            StringToken st = NM_StringCheck->addString(first_kmer);
            std::cout << "SP: " << working_str<<" : "<<first_kmer<<"("<<st<<")"<<std::endl;
            
            CrisprNode * first_node = new CrisprNode(st);
            
            first_node->addPartner(prev_node, true);
            NM_Nodes.push_back(first_node);
        }
    }
    else
    {
        // start with a spacer
        
        // for the first spacer, cut only the second kmer ( if starting with a spacer )
        // as there is not a direct repeat to anchor the other side
        if(RH->getFirstSpacer(&working_str))
        {
            // make sure that we can cut a kmer and add the node
            if (CRASS_DEF_NODE_KMER_SIZE < working_str.length()) 
            {
                std::string second_kmer = working_str.substr(working_str.length() - CRASS_DEF_NODE_KMER_SIZE, CRASS_DEF_NODE_KMER_SIZE);
                StringToken st = NM_StringCheck->addString(second_kmer);
                std::cout << "SP: " << working_str<<" : "<<second_kmer<<"("<<st<<")"<<std::endl;
                
                CrisprNode * first_node = new CrisprNode(st);
                
                NM_Nodes.push_back(first_node);
                
                prev_node = first_node;
            }

            size_of_RH_startStops -= 2;
        }
        else
        {
            logError("could't get the first spacer");
        }
        
        while(size_of_RH_startStops > 0)
        {

            if(RH->getNextSpacer(&working_str))
            {
                addCrisprNodes(&prev_node, working_str);
                size_of_RH_startStops -= 2;

            }
            else
            {
                break;
            }
        }
        
        // do the last spacer if required
        // TODO  -1?
        if (RH->back() <= RH->seqLength() - CRASS_DEF_NODE_KMER_SIZE) 
        {
            RH->getNextSpacer(&working_str);
            
            std::string first_kmer = working_str.substr(0, CRASS_DEF_NODE_KMER_SIZE);
            StringToken st = NM_StringCheck->addString(first_kmer);
            std::cout << "SP: " << working_str<<" : "<<first_kmer<<"("<<st<<")"<<std::endl;
            
            CrisprNode * first_node = new CrisprNode(st);
            first_node->addPartner(prev_node, true);

            NM_Nodes.push_back(first_node);
        }
    }    
}


//----
// Private function called from splitReadHolder to cut the kmers and make the nodes
//
void NodeManager::addCrisprNodes(CrisprNode ** prevNode, std::string& workingString)
{
    // now cut kmers on either side and add the pair into the node manager 
    std::string first_kmer = workingString.substr(0,CRASS_DEF_NODE_KMER_SIZE);
    std::string second_kmer = workingString.substr(workingString.length() - CRASS_DEF_NODE_KMER_SIZE, CRASS_DEF_NODE_KMER_SIZE );
    
    
    StringToken st1 = NM_StringCheck->addString(first_kmer);
    StringToken st2 = NM_StringCheck->addString(second_kmer);
    
    std::cout << "SP: " << workingString<<" : "<<first_kmer<<"("<<st1<<") : "<<second_kmer<<"("<<st2<<")"<<std::endl;

    CrisprNode * first_kmer_node = new CrisprNode(st1);
    CrisprNode * second_kmer_node = new CrisprNode(st2);
    
    
    
    // the first kmers pair is the previous node which lay before it therefore bool is true
    // make sure prevNode is not NULL
    if (*prevNode != NULL) 
    {
        first_kmer_node->addPartner(*prevNode, true);
    }
    
    second_kmer_node->addPartner(first_kmer_node, true);
    
    NM_Nodes.push_back(first_kmer_node);
    NM_Nodes.push_back(second_kmer_node);

    StringToken sp_str_token = NM_StringCheck->addString(workingString);
    
    // create a spacer instance from these two nodes
    SpacerInstance * curr_spacer = new SpacerInstance(sp_str_token, first_kmer_node, second_kmer_node); 
    
    NM_Spacers.push_back(curr_spacer);
    
    *prevNode = second_kmer_node;
}


//----
// Private function called from addCrisprNodes: adds in a spacer instance
//
void NodeManager::addSpacerInstance(SpacerInstance * newSpacer)
{
    NM_Spacers.push_back(newSpacer);
}






















