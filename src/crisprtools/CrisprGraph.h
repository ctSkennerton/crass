/*
 *  CrisprGraph.h is part of the crisprtools project
 *  
 *  Created by Connor Skennerton on 7/12/11.
 *  Copyright 2011 Connor Skennerton. All rights reserved. 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
 *                     A B R A K A D A B R A
 *                      A B R A K A D A B R
 *                       A B R A K A D A B
 *                        A B R A K A D A       	
 *                         A B R A K A D
 *                          A B R A K A
 *                           A B R A K
 *                            A B R A
 *                             A B R
 *                              A B
 *                               A
 */

#ifndef crisprtools_CrisprGraph_h
#define crisprtools_CrisprGraph_h

#include <graphviz/gvc.h>
#include <map>

namespace crispr {
    
    // basic wrapper around the graphviz library
    class graph {
      
        Agraph_t * G_graph;
        
    public:
        graph(char * name )
        {
            G_graph = agopen(name, AGDIGRAPH );
        }
        
        ~graph()
        {
            agclose(G_graph);
        }
        
        inline Agraph_t * getGraph(void){return G_graph;}
        
        
        inline Agnode_t * addNode(char * name)
        {
            Agnode_t * n = agnode(G_graph, name);
            return n;
        }       
        
        inline Agedge_t * addEdge(Agnode_t * tail, Agnode_t * head)
        {
            return agedge(G_graph, tail, head);
        } 
        
        inline Agnode_t * nodeWithName(char * name)
        {
            return agfindnode(G_graph, name);
        }      
        
        inline Agedge_t * findEdge(Agnode_t * tail, Agnode_t * head)
        {
            return agfindedge(G_graph, tail, head);
        }
        
        inline int setNodeAttribute(Agnode_t * node ,char * attrName, char * attrValue)
        {
            return agsafeset(node, attrName, attrValue, "");
        }  
        inline int setGraphAttribute(char * attrName, char * attrValue)
        {
            return agsafeset(G_graph, attrName, attrValue, "");
        }
        inline int setEdgeAttribute(Agedge_t * edge ,char * attrName, char * attrValue)
        {
            return agsafeset(edge, attrName, attrValue, "");
        }
        inline int setNodeAttribute(char * nodeName ,char * attrName, char * attrValue)
        {
            Agnode_t * node = agfindnode(G_graph, nodeName);
            if (node != NULL) {
                return setNodeAttribute(node, attrName, attrValue);
            } else {
                return 1;
            }
        }

    };
}


#endif
