// DrawTool.cpp
//
// Copyright (C) 2011 - Connor Skennerton
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "DrawTool.h"
#include <libcrispr/Exception.h>
#include "CrisprGraph.h"
#include "Utils.h"
#include "config.h"
#include <libcrispr/StlExt.h>
#include <string.h>
#include <sys/stat.h>
#include <graphviz/gvc.h>
#include <getopt.h>

DrawTool::~DrawTool()
{
    graphVector::iterator iter = DT_Graphs.begin();
    while (iter != DT_Graphs.end()) {
        if ((*iter) != NULL) {
            delete *iter;
            *iter = NULL;
        }
        iter++;
    }
    gvFreeContext(DT_Gvc);

}

int DrawTool::processOptions (int argc, char ** argv)
{
	try {
        int c;
        int index;
        static struct option long_options [] = {       
            {"help", no_argument, NULL, 'h'},
            {"outfile", required_argument, NULL, 'o'},
            {"colour",required_argument,NULL,'c'},
            {"bins", required_argument, NULL, 'b'},
            {"format", required_argument, NULL, 'f'},
            {"algorithm", required_argument, NULL, 'a'},
            {"groups", required_argument, NULL, 'g'},
            {0,0,0,0}
        };
        
        bool algo = false, outformat = false;
        while((c = getopt_long(argc, argv, "hg:c:a:f:o:b:", long_options, &index)) != -1)
        {
            switch(c)
            {
                case 'h':
                {
                    drawUsage ();
                    exit(1);
                    break;
                }
                case 'g':
                {
                    generateGroupsFromString (optarg);
                    break;
                }
                    
                case 'o':
                {
                    DT_OutputFile = optarg;
                    if (DT_OutputFile[DT_OutputFile.length() - 1] != '/')
                    {
                        DT_OutputFile += '/';
                    }
                    
                    // check if our output folder exists
                    struct stat file_stats;
                    if (0 != stat(DT_OutputFile.c_str(),&file_stats)) 
                    {
                        recursiveMkdir(DT_OutputFile);
                    }
                    break;
                }
                case 'c':
                {
                    if (!strcmp(optarg, "red-blue") ) {
                        DT_ColourType = RED_BLUE;
                    } else if (!strcmp(optarg, "red-blue-green")) {
                        DT_ColourType = RED_BLUE_GREEN;
                    } else if (!strcmp(optarg, "blue-red")) {
                        DT_ColourType = BLUE_RED;
                    } else if (!strcmp(optarg, "green-blue-red")) {
                        DT_ColourType = GREEN_BLUE_RED;
                    } else {
                        throw crispr::input_exception("Not a known color type");
                    }
                    break;
                } 
                case 'f':
                {
                    DT_OutputFormat = optarg;
                    outformat = true;
                    break;
                }
                case 'a':
                {
                    if (!strcmp(optarg, "dot") || 
                        !strcmp(optarg, "neato") || 
                        !strcmp(optarg, "fdp") || 
                        !strcmp(optarg, "sfdp") ||
                        !strcmp(optarg, "twopi") || 
                        !strcmp(optarg, "circo")) {
                        DT_RenderingAlgorithm = optarg;       
                    } else {
                        throw crispr::input_exception("Not a known Graphviz rendering algorithm");
                    }
                    algo = true;
                    break;
                }
                case 'b':
                {
                    int i;
                    if (from_string<int>(i, optarg, std::dec)) {
                        if (i > 0) {
                            DT_Bins = i;
                        } else {
                            throw crispr::input_exception("The number of bins of colour must be greater than 0");
                        }
                    } else {
                        throw crispr::runtime_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,"could not convert string to int");
                    }
                    break;
                }
                default:
                {
                    drawUsage();
                    exit(1);
                    break;
                }
            }
        }
        if (!(algo & outformat) ) {
            throw crispr::input_exception("You must specify both -a and -f on the command line");
        }

    } catch (crispr::input_exception& e) {
        std::cerr<<e.what()<<std::endl;
        drawUsage();
        exit(1);
    } catch (crispr::runtime_exception& e) {
        std::cerr<<e.what()<<std::endl;
        exit(1);
    }
    return optind;
}

void DrawTool::generateGroupsFromString ( std::string str)
{
	std::set<std::string> vec;
	split ( str, vec, ",");
	DT_Groups = vec;
	DT_Subset = true;
}

int DrawTool::processInputFile(const char * inputFile)
{
    try {
        crispr::xml::parser xml_parser;
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        
        if( !root_elem ) throw(crispr::xml_exception(__FILE__, 
                                                     __LINE__, 
                                                     __PRETTY_FUNCTION__, 
                                                     "empty XML document" ));

        // get the children
        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); 
             currentElement != NULL; 
             currentElement = currentElement->getNextElementSibling()) {
                
            // is this a group element
            if (xercesc::XMLString::equals(currentElement->getTagName(), xml_parser.tag_Group())) {
                char * c_group_id = tc(currentElement->getAttribute(xml_parser.attr_Gid()));
                std::string group_id = c_group_id;
                if (DT_Subset) {
                    
                    // we only want some of the groups look at DT_Groups
                    if (DT_Groups.find(group_id.substr(1)) != DT_Groups.end() ) {
                        parseGroup(currentElement, xml_parser);
                        
                    }
                } else {
                    parseGroup(currentElement, xml_parser);   
                }
                xr(&c_group_id);
            }
            
        }
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    } catch (crispr::runtime_exception& e) {
        std::cerr<<e.what()<<std::endl;
    }
    
    return 0;
}
void DrawTool::parseGroup(xercesc::DOMElement * parentNode, crispr::xml::parser& xmlParser)
{

    
    // create a new graph object
    char * c_gid = tc(parentNode->getAttribute(xmlParser.attr_Gid()));
    crispr::graph * current_graph = new crispr::graph(c_gid);
    xr(&c_gid);
    // change the max and min coverages back to their original values
    resetInitialLimits();
    
    DT_Graphs.push_back(current_graph);
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Data())) {
            parseData(currentElement, xmlParser, current_graph);
            setColours();
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Assembly())) {
            parseAssembly(currentElement, xmlParser, current_graph);
        }
        
    }
    
    char * c_file_prefix = tc(parentNode->getAttribute(xmlParser.attr_Gid()));
    std::string file_prefix = c_file_prefix;
    std::string file_name = file_prefix + "." + DT_OutputFormat;
    xr(&c_file_prefix);
    char * file_name_c = strdup(file_name.c_str());

    layoutGraph(current_graph->getGraph(), DT_RenderingAlgorithm);
    renderGraphToFile(current_graph->getGraph(), DT_OutputFormat, file_name_c);
    freeLayout(current_graph->getGraph());
    // free the duplicated string
    try {
        delete file_name_c;
    } catch (std::exception& e) {
        std::cerr<<e.what()<<std::endl;
    }
}

void DrawTool::parseData(xercesc::DOMElement * parentNode, 
                         crispr::xml::parser& xmlParser, 
                         crispr::graph * currentGraph)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        /*if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Drs())) {
                // change the direct repeats
                parseDrs(currentElement, xmlParser);
        } else*/
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacers())) {
                // change the spacers
                parseSpacers(currentElement, xmlParser, currentGraph);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flankers())) {
                // change the flankers
                parseFlankers(currentElement, xmlParser, currentGraph);
        }
    }
}


//void DrawTool::parseDrs(xercesc::DOMElement * parentNode, crispr::xml::parser& xmlParser)
//{
//    xercesc::DOMNodeList * children = parentNode->getChildNodes();
//    const  XMLSize_t nodeCount = children->getLength();
//    
//    // For all nodes, children of "root" in the XML tree.
//    for( XMLSize_t xx = 0; xx < nodeCount; ++xx ) {
//        xercesc::DOMNode * currentNode = children->item(xx);
//        if( currentNode->getNodeType() &&  currentNode->getNodeType() == xercesc::DOMNode::ELEMENT_NODE ) {
//            // Found node which is an Element. Re-cast node as element
//            xercesc::DOMElement* currentElement = dynamic_cast< xercesc::DOMElement* >( currentNode );
//            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getDr())) {
//
//            }
//        }
//    }
//}
void DrawTool::parseSpacers(xercesc::DOMElement * parentNode, 
                            crispr::xml::parser& xmlParser, 
                            crispr::graph * currentGraph)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacer())) {
            char * node_name = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            
            Agnode_t * current_graphviz_node = currentGraph->addNode(node_name);
            
            char * shape = strdup("shape");
            char * circle = strdup("circle");
            currentGraph->setNodeAttribute(current_graphviz_node, shape, circle);
            delete shape;
            delete circle;
            
            
            char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            std::string spid = c_spid;

            if (currentElement->hasAttribute(xmlParser.attr_Cov())) {
                char * c_cov = tc(currentElement->getAttribute(xmlParser.attr_Cov()));
                double current_cov;
                if (from_string<double>(current_cov, c_cov, std::dec)) {
                    recalculateLimits(current_cov);
                    DT_SpacerCoverage[spid] = std::pair<bool, double>(true,current_cov);
                } else {
                    throw crispr::runtime_exception(__FILE__, 
                                                    __LINE__, 
                                                    __PRETTY_FUNCTION__,
                                                    "Unable to convert serialized coverage");
                }
                
                xr(&c_cov);
                
            } else {
                DT_SpacerCoverage[spid] = std::pair<bool, double>(false,0);
            }
            
            xr(&c_spid);
            xr(&node_name);
        }
        
    }
}

void DrawTool::parseFlankers(xercesc::DOMElement * parentNode, 
                             crispr::xml::parser& xmlParser, 
                             crispr::graph * currentGraph)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flanker())) {
            char * c_flid = tc(currentElement->getAttribute(xmlParser.attr_Flid()));
            Agnode_t * current_graphviz_node = currentGraph->addNode(c_flid);
            
            char * shape = strdup("shape");
            char * shape_val = strdup("diamond");
            
            currentGraph->setNodeAttribute(current_graphviz_node, shape, shape_val);
            
            delete shape;
            delete shape_val;
            xr(&c_flid);
        }
    }
}

void DrawTool::parseAssembly(xercesc::DOMElement * parentNode, 
                             crispr::xml::parser& xmlParser, 
                             crispr::graph * currentGraph)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Contig())) {
            char * c_contig_id = tc(currentElement->getAttribute(xmlParser.attr_Cid()));
            std::string contig_id = c_contig_id;
            parseContig(currentElement, xmlParser,currentGraph, contig_id);
            xr(&c_contig_id);
        }
            
        
    }
    
}

void DrawTool::parseContig(xercesc::DOMElement * parentNode, 
                           crispr::xml::parser& xmlParser, 
                           crispr::graph * currentGraph, 
                           std::string& contigId)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Cspacer())) {
            // get the node
            char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            Agnode_t * current_graphviz_node = currentGraph->addNode(c_spid);

            // find the coverage for this guy if it was set

            std::string cov_spid = c_spid;
			xr(&c_spid);
            if (DT_SpacerCoverage[cov_spid].first) {
                std::string color = DT_Rainbow.getColour(DT_SpacerCoverage[cov_spid].second);

                // fix things up for Graphviz
                char * color_for_graphviz = strdup(('#' + color).c_str());

                // add in the colour attributes
                char * style = strdup("style");
                char * filled = strdup("filled");
                char * fillcolour = strdup("fillcolor");
                
                currentGraph->setNodeAttribute(current_graphviz_node, style, filled);
                currentGraph->setNodeAttribute(current_graphviz_node, fillcolour, color_for_graphviz );
                
                delete style;
                delete filled;
                delete fillcolour;
                delete color_for_graphviz;
            }

            parseCSpacer(currentElement, xmlParser,currentGraph,current_graphviz_node,contigId);
        }
    }
}

void DrawTool::parseCSpacer(xercesc::DOMElement * parentNode, 
                            crispr::xml::parser& xmlParser, 
                            crispr::graph * currentGraph, 
                            Agnode_t * currentGraphvizNode, 
                            std::string& contigId)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

                    
       /* if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getBspacers())) {
                parseLinkSpacers(currentElement, xmlParser,currentGraph,currentGraphvizNode, REVERSE,contigId);
        } else*/ if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Fspacers())) {
                parseLinkSpacers(currentElement, xmlParser,currentGraph,currentGraphvizNode, FORWARD,contigId);
        /*} else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getBflankers())) {
                parseLinkFlankers(currentElement, xmlParser,currentGraph,currentGraphvizNode, REVERSE,contigId);
        */} else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Fflankers())) {
                parseLinkFlankers(currentElement, xmlParser,currentGraph,currentGraphvizNode, FORWARD,contigId);
        }
        
    }
}

void DrawTool::parseLinkSpacers(xercesc::DOMElement * parentNode, 
                                crispr::xml::parser& xmlParser, 
                                crispr::graph * currentGraph, 
                                Agnode_t * currentGraphvizNode, 
                                EDGE_DIRECTION edgeDirection, 
                                std::string& contigId)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
            
        char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
        Agnode_t * edge_node = currentGraph->addNode(c_spid);
        xr(&c_spid);
        if (edgeDirection == FORWARD) {
            currentGraph->addEdge(currentGraphvizNode, edge_node);
        } else {
            currentGraph->addEdge(edge_node, currentGraphvizNode);
        }
    }
}

void DrawTool::parseLinkFlankers(xercesc::DOMElement * parentNode, 
                                 crispr::xml::parser& xmlParser, 
                                 crispr::graph * currentGraph, 
                                 Agnode_t * currentGraphvizNode, 
                                 EDGE_DIRECTION edgeDirection, 
                                 std::string& contigId)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        char * c_flid = tc( currentElement->getAttribute(xmlParser.attr_Flid()));

        Agnode_t * edge_node = currentGraph->addNode(c_flid);
        xr(&c_flid);
        if (edgeDirection == FORWARD) {
            currentGraph->addEdge(currentGraphvizNode, edge_node);
        } else {
            currentGraph->addEdge(edge_node, currentGraphvizNode);
        }        
    }
}


int drawMain (int argc, char ** argv)
{
    try {
		DrawTool dt;
		int opt_index = dt.processOptions (argc, argv);
		if (opt_index >= argc) {
			throw crispr::input_exception("No input file provided" );
            
        } else {
			// get cracking and process that file
			return dt.processInputFile(argv[opt_index]);
		}
	} catch(crispr::input_exception& re) {
        std::cerr<<re.what()<<std::endl;
        drawUsage();
        return 1;
    } catch(crispr::exception& ce ) {
		std::cerr<<ce.what()<<std::endl;
		return 1;
	}
}

void drawUsage(void)
{
    std::cout<<PACKAGE_NAME<<" draw [-ghyo] -a ALGORITHM -f FORMAT file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h					print this handy help message"<<std::endl;
    std::cout<<"-o DIR              output file directory  [default: .]" <<std::endl; 
	std::cout<<"-g INT[,n]          a comma separated list of group IDs that you would like to extract data from."<<std::endl;
	std::cout<<"					Note that only the group number is needed, do not use prefixes like 'Group' or 'G', which"<<std::endl;
	std::cout<<"					are sometimes used in file names or in a .crispr file"<<std::endl;
	std::cout<<"-a STRING           The Graphviz layout algorithm to use [default: dot ]"<<std::endl;
    std::cout<<"-f STRING           The output format for the image, equivelent to the -T parameter of Graphviz executables [default: eps]"<<std::endl;
    std::cout<<"-b INT              Number of colour bins"<<std::endl;
    std::cout<<"-c COLOUR           The colour scale to use for coverage information.  The available choices are:"<<std::endl;
    std::cout<<"                        red-blue"<<std::endl;
    std::cout<<"                        blue-red"<<std::endl;
    std::cout<<"                        red-blue-green"<<std::endl;
    std::cout<<"                        green-blue-red"<<std::endl;
}
