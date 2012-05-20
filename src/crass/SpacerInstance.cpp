// File: SpacerInstance.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// Implementation of SpacerInstance functions
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
#include <string>

// local includes
#include "SpacerInstance.h"
#include "CrisprNode.h"
#include "StringCheck.h"
#include "LoggerSimp.h"

SpacerInstance::SpacerInstance(StringToken spacerID) 
{
    //-----
    // constructor
    //
    SI_SpacerSeqID = spacerID;
    SI_LeadingNode = NULL;
    SI_LastNode = NULL;
    SI_InstanceCount = 0;
    SI_ContigID = 0;
    SI_Attached = false;
    SI_isFlanker = false;
}

SpacerInstance::SpacerInstance(StringToken spacerID, CrisprNode * leadingNode, CrisprNode * lastNode)
{
    SI_SpacerSeqID = spacerID;
    SI_LeadingNode = leadingNode;
    SI_LastNode = lastNode;
    SI_InstanceCount  = 1;
    SI_ContigID = 0;
    SI_Attached = false;
    SI_isFlanker = false;
}

void SpacerInstance::clearEdge(void)
{
	//-----
	// Clear all edges
	//
    SpacerEdgeVector_Iterator iter = SI_SpacerEdges.begin();
    while (iter != SI_SpacerEdges.end()) 
    {
        if (*iter != NULL) 
        {
            delete *iter;
            *iter = NULL;
        }
        iter++;
    }
}

bool SpacerInstance::isFur(void)
{
	//-----
	// Check to see if the spacer is a cap joined onto a non-cap
	//
	SpacerEdgeVector_Iterator edge_iter = SI_SpacerEdges.begin();
	
	// zero rank spacers are viable by default
	if(1 != getSpacerRank())
		return false;
			
	while(edge_iter != SI_SpacerEdges.end())
	{
		if((*edge_iter)->edge->getSpacerRank() > 2)
		{
			return true;
		}
		edge_iter++;
	}
	return false;
}


bool SpacerInstance::isViable(void)
{
	//-----
	// Check to see if the spacer has forward AND reverse edges
	// return true or false accordingly
	//
	SpacerEdgeVector_Iterator edge_iter = SI_SpacerEdges.begin();
	
	// zero rank spacers are viable by default

	if(getSpacerRank() < 2) {
    
		return true;
    }
			
	bool has_forward = false;
	bool has_reverse = false;
	while(edge_iter != SI_SpacerEdges.end())
	{
        if((*edge_iter)->d == REVERSE)
			has_reverse = true;
		else
			has_forward = true;
		if(has_reverse && has_forward)
			return true;
		edge_iter++;
	}
	return false;
}

SpacerEdgeVector_Iterator SpacerInstance::find(SpacerInstance * si)
{
    // check to see if the spacer instance in in the edges
    SpacerEdgeVector_Iterator iter;
    for (iter = SI_SpacerEdges.begin(); iter != SI_SpacerEdges.end(); ++iter) {
        if ((*iter)->edge == si) {
            return iter;
        }
    }
    return SI_SpacerEdges.end();
}

void SpacerInstance::detachFromSpacerGraph(void)
{
	//-----
	// remove this spacer from the graph
	//
	if(0 == getSpacerRank())
		return;
	SpacerEdgeVector_Iterator edge_iter = SI_SpacerEdges.begin();
	while(edge_iter != SI_SpacerEdges.end())
	{
		// delete the return edge
		//std::cout << "c: " << this << " : " << (*edge_iter)->edge << std::endl;
		if(!(*edge_iter)->edge->detachSpecificSpacer(this))
		{
			return;
		}
        
		// free the memory!
		if (*edge_iter != NULL) 
        {
            delete *edge_iter;
            *edge_iter = NULL;
        }
		edge_iter++;
	}
	
	// no edges left here
	SI_SpacerEdges.clear();
}

bool SpacerInstance::detachSpecificSpacer(SpacerInstance * target)
{
	//-----
	// remove this spacer from the graph
	//
	if(0 == getSpacerRank())
	{
		logError("Trying to remove edge from zero rank spacer");
		return false;
	}
	
	SpacerEdgeVector_Iterator edge_iter = SI_SpacerEdges.begin();

	while(edge_iter != SI_SpacerEdges.end())
	{
		// we need to find the target edge
		if((*edge_iter)->edge == target)
		{
			//std::cout << "rev: " << target << " : " << this << std::endl;
			// free the memory
	        if (*edge_iter != NULL) 
	        {
	            delete *edge_iter;
	            *edge_iter = NULL;
	        }
	        
	        // remove from the vector
	        SI_SpacerEdges.erase(edge_iter);
	        return true;
		}
		
		edge_iter++;
	}
	
	logError("Could not find target: " << target);
	return false;
}

void SpacerInstance::printContents(void)
{
	//-----
	// print the contents of all the contents
	//
	
	std::cout << "-------------------------------\n" << this << std::endl;
	std::cout << "ST: " << SI_SpacerSeqID << " LEADER: " << SI_LeadingNode->getID() << " LAST: " << SI_LastNode->getID() << std::endl;
	std::cout << "IC: " << SI_InstanceCount << " ATT? " << SI_Attached << " CID: " << SI_ContigID << std::endl;
	SpacerEdgeVector_Iterator edge_iter = SI_SpacerEdges.begin();
	while(edge_iter != SI_SpacerEdges.end())
	{
		std::cout << " --> " << (*edge_iter)->edge << " : " << (*edge_iter)->d << std::endl;
		edge_iter++;
	}
}



