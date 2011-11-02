/*
 *  GraphDrawingDefines.h is part of the crass project
 *  
 *  Created by Connor Skennerton on 6/09/11.
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

#ifndef crass_GraphDrawingDefines_h
#define crass_GraphDrawingDefines_h

// --------------------------------------------------------------------
// GRAPH DRAWING DEFAULTS
// --------------------------------------------------------------------
#define CRASS_DEF_GV_DEFULT_COLOUR              "grey"              // colour when no coverage is available 
#define CRASS_DEF_GV_EDGE_LENGTH                2                   // length of arrows for edges in images
#define CRASS_DEF_GV_NODE_SHAPE                 "circle"            // shape of nodes
#define CRASS_DEF_GV_NODE_PREFIX                "node_"             // prefix used when naming nodes for the image
#define CRASS_DEF_GV_NODE_FILL                  "filled"            // should our nodes be filled with colour
#define CRASS_DEF_DEFAULT_GRAPH_TYPE            "digraph"           // the default type of graphviz graph to generate

#define CRASS_DEF_DEFAULT_SUB_GRAPH_TYPE 		"subgraph"           // the default type of graphviz graph to generate
#define CRASS_DEF_GV_SPA_EDGE_LENGTH                2                   // length of arrows for edges in images
#define CRASS_DEF_GV_SPA_SHAPE                 "rectangle"            // shape of nodes
#define CRASS_DEF_GV_SPA_PREFIX                 "sp_"               // prefix used when naming spacers for the image

// --------------------------------------------------------------------
// GRAPH DRAWING MACROS
// --------------------------------------------------------------------
#define gvEdge(stream,id1,id2){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" -> "<<CRASS_DEF_GV_NODE_PREFIX<<id2<< " [ len="<<CRASS_DEF_GV_EDGE_LENGTH<<" ];"<<std::endl;\
}

#define gvJumpingEdge(stream,id1,id2){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" -> "<<CRASS_DEF_GV_NODE_PREFIX<<id2<< " [ len="<<CRASS_DEF_GV_EDGE_LENGTH<<", style=dashed ];"<<std::endl;\
}

#define gvCustomEdge(stream,id1,id2, cUSTOMsTRING){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" -> "<<CRASS_DEF_GV_NODE_PREFIX<<id2<< " [ "<<cUSTOMsTRING<<" ];"<<std::endl;\
}

#define gvNodeF(stream,id1,color){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" [ color = \"#"<<color<<"\", fillcolor=\"#"<<color<<"\", style= "<<CRASS_DEF_GV_NODE_FILL<<\
                                                ", shape="<< CRASS_DEF_GV_NODE_SHAPE<<"];"<<std::endl;\
}

#define gvNodeB(stream,id1,color){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" [ color = \"#"<<color<<"\", fillcolor=\"#"<<color<<"\", style= "<<CRASS_DEF_GV_NODE_FILL<<\
                                                ", shape="<< CRASS_DEF_GV_NODE_SHAPE<<"];"<<std::endl;\
}

#define gvCustomNode(stream,id1,cUSTOMsTRING){\
    stream<<CRASS_DEF_GV_NODE_PREFIX<<id1<<" [ "<<cUSTOMsTRING<<" ];"<<std::endl;\
}

#define gvGraphHeader(stream,gRAPHtITLE){\
    stream<<CRASS_DEF_DEFAULT_GRAPH_TYPE<<" "<<gRAPHtITLE<<" {"<<std::endl;\
}

#define gvGraphFooter(stream){ \
    stream<<std::endl<<"}"<<std::endl;\
}

#define gvSpacer(stream,id1,color){\
    stream<<"\t\t" << id1<<" [ color = \"#"<<color<<"\", fillcolor=\"#"<<color<<"\", style= "<<CRASS_DEF_GV_NODE_FILL<<\
                                                ", shape="<< CRASS_DEF_GV_SPA_SHAPE<<"];"<<std::endl;\
}

#define gvSpEdge(stream,id1,id2){\
    stream<<"\t\t"<<id1<<" -> "<<id2<< " [ len="<<CRASS_DEF_GV_SPA_EDGE_LENGTH<<" ];"<<std::endl;\
}

#define gvSubGraphHeader(stream,gRAPHtITLE,cLUSTERnUMBER){\
    stream<<"\t"<<CRASS_DEF_DEFAULT_SUB_GRAPH_TYPE<<" cluster_"<<cLUSTERnUMBER<<"\t{\n\t\tstyle=filled;\n\t\tcolor=lightgrey;\n\t\tlabel = \""<<gRAPHtITLE <<"\";\n"<<std::endl;\
}

#define gvSubGraphFooter(stream){ \
    stream<<std::endl<<"\t}"<<std::endl;\
}

#define gvKeyHeader(stream, cLUSTERnUMBER){ \
    stream<<"\t"<<CRASS_DEF_DEFAULT_SUB_GRAPH_TYPE<<" cluster_"<<cLUSTERnUMBER<<"\t{\n\t\t\"KEY\" [ fillcolor = \"white\" shape = \"record\" label =<<table border=\"0\" cellborder=\"0\" cellpadding=\"0\" bgcolor=\"white\"><tr><td>KEY</td></tr>";\
}

#define gvKeyEntry(stream, lABEL, cOLOUR){ \
    stream<<"<tr><td bgcolor=\"#" << cOLOUR << "\" align=\"center\" colspan=\"2\"><font color=\"white\">"<< lABEL <<"</font></td></tr>";\
}

#define gvKeyFooter(stream){ \
    stream<<"</table>> ];\n\t}"<<std::endl;\
}

#endif
