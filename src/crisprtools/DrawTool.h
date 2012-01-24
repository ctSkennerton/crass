/*
 * DrawTool.h
 *
 * Copyright (C) 2011 - Connor Skennerton
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "XML.h"
#include "CrisprGraph.h"
#include "Rainbow.h"
#include <graphviz/gvc.h>
#include <set>
#include <string>
#include <map>

#ifndef DRAWTOOL_H
#define DRAWTOOLS_H 
typedef std::vector<crispr::graph * > graphVector;
class DrawTool {
    
    enum EDGE_DIRECTION {
        FORWARD = 0,
        REVERSE = 1
    };
    
    GVC_t * DT_Gvc;
    std::string DT_OutputFile;
    char * DT_RenderingAlgorithm;
    char * DT_OutputFormat;
    std::set<std::string> DT_Groups;
    std::map<std::string, std::pair<bool, double> > DT_SpacerCoverage;
    bool DT_Subset;
    graphVector DT_Graphs;
    Rainbow DT_Rainbow;
    RB_TYPE DT_ColourType;
    int DT_Bins;
    double DT_UpperLimit;
    double DT_LowerLimit;
    
    
    void resetInitialLimits(void) 
    {        
        DT_LowerLimit = 10000000;
        DT_UpperLimit = 0;
    }
    
    void recalculateLimits(double d) 
    {
        if (d > DT_UpperLimit) 
        {
            DT_UpperLimit = d;
        }
        if(d < DT_LowerLimit)
        {
            DT_LowerLimit = d;
        }
    }
    
    void setColours(void)
    {
        DT_Rainbow.setType(DT_ColourType);
        if (DT_Bins > 0) 
        {
            DT_Rainbow.setLimits(DT_LowerLimit, DT_UpperLimit, DT_Bins);
        } 
        else 
        {
            DT_Rainbow.setLimits(DT_LowerLimit, DT_UpperLimit);
        }
    }

public:
    DrawTool()
    {
        DT_Gvc = gvContext();
        DT_Subset = false;
        DT_RenderingAlgorithm = "dot";
        DT_OutputFormat = "eps";
        DT_LowerLimit = 10000000;
        DT_UpperLimit = 0;
        DT_ColourType = RED_BLUE;
        DT_Bins = -1;
    }
    
    ~DrawTool();
    
    
    inline GVC_t * getContext(void){return DT_Gvc;}
    
    void layoutGraph(Agraph_t * g, char * a){gvLayout(DT_Gvc, g, a);}
    void renderGraphToFile(Agraph_t * g, char * t, char * f){gvRenderFilename(DT_Gvc, g, t, f);}
    void freeLayout(Agraph_t * g){gvFreeLayout(DT_Gvc, g);}
    
    int processOptions(int argc, char ** argv);
    void generateGroupsFromString ( std::string str);
    int processInputFile(const char * inputFile);
    void parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * current_graph);
    void parseDrs(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * current_graph);
    void parseSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * current_graph);
    void parseFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * current_graph);
    
    void parseAssembly(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * currentGraph);
    void parseContig(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * currentGraph, std::string& contigId);
    void parseCSpacer(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * currentGraph, Agnode_t * currentGraphvizNode, std::string& contigId);
    void parseLinkSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * currentGraph, Agnode_t * currentGraphvizNode, EDGE_DIRECTION edgeDirection, std::string& contigId);
    void parseLinkFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, crispr::graph * currentGraph, Agnode_t * currentGraphvizNode, EDGE_DIRECTION edgeDirection, std::string& contigId);

};

int drawMain(int argc, char ** argv);
void drawUsage(void);
#endif
