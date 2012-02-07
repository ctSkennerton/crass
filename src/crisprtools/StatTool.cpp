// StatTool.cpp
//
// Copyright (C) 2011 - Connor Skennerton
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "StatTool.h"
#include "StlExt.h"
#include "config.h"
#include "Exception.h"
#include <iostream>
#include <getopt.h>
#include <cstring>

StatTool::~StatTool()
{
    std::vector<StatManager *>::iterator iter = begin();
    while (iter != end()) {
        delete *iter;
        iter++;
    }
}

void StatTool::generateGroupsFromString ( std::string str)
{
	std::set<std::string> vec;
	split ( str, vec, ",");
	ST_Groups = vec;
	ST_Subset = true;
}
int StatTool::processOptions (int argc, char ** argv)
{
	int c;
	while((c = getopt(argc, argv, "ahHg:ps:to:")) != -1)
	{
        switch(c)
		{
			case 'a':
            {
                ST_AggregateStats = true;
                break;
            }
            case 'h':
			{
				statUsage();
				exit(1);
				break;
			}
            case 'o':
            {
                ST_OutputFileName = optarg;
                break;
            }
            case 'p':
            {
                ST_Pretty = true;
                break;
            }
            case 'g':
            {
                generateGroupsFromString(optarg);
                break;
            }
            case 's':
            {
                ST_Separator = optarg;
                break;
            }
            case 't':
            {
                ST_Tabular = true;
                break;
            }
            case 'H':
            {
                ST_WithHeader = true;
                break;
            }
            default:
            {
                statUsage();
                exit(1);
                break;
            }
		}
        try {
            if(!ST_Pretty && !ST_Tabular && !ST_AggregateStats) {
                throw crispr::input_exception("You must specify one of -t -a -p");
            }
            if(ST_Pretty && ST_Tabular) {
                throw crispr::input_exception("-p and -t cannot both be specified");
            }
        } catch (crispr::input_exception& e) {
            std::cerr<< e.what()<<std::endl;
            statUsage();
            exit(1);
        }
	}
	return optind;
}

int StatTool::processInputFile(const char * inputFile)
{
    try {
        crispr::XML xml_parser;
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        int num_groups_to_process = static_cast<int>(ST_Groups.size());
        //std::cout<<num_groups_to_process<<std::endl;
        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {

            if (ST_Subset && num_groups_to_process == 0) {
                break;
            }
            // is this a group element
            if (xercesc::XMLString::equals(currentElement->getTagName(), xml_parser.getGroup())) {
                char * c_gid = tc(currentElement->getAttribute(xml_parser.getGid()));
                std::string group_id = c_gid;
                if (ST_Subset) {
                    // we only want some of the groups look at DT_Groups
                    if (ST_Groups.find(group_id.substr(1)) != ST_Groups.end() ) {
                        parseGroup(currentElement, xml_parser);

                        // decrease the number of groups left
                        // if we are only using a subset
                        if(ST_Subset) num_groups_to_process--;
                    }
                } else {
                    parseGroup(currentElement, xml_parser);   
                }
                xr(&c_gid);
            }
        }
        AStats agregate_stats;
        agregate_stats.total_groups = 0;
        agregate_stats.total_spacers = 0;
        agregate_stats.total_dr = 0;
        agregate_stats.total_flanker = 0;
        agregate_stats.total_spacer_length = 0;
        agregate_stats.total_spacer_cov = 0;
        agregate_stats.total_dr_length = 0;
        agregate_stats.total_flanker_length = 0;
        // go through each of the groups and print out a pretty picture
        std::vector<StatManager *>::iterator iter = this->begin();
        while (iter != this->end()) {
            if (ST_Pretty) {
                prettyPrint(*iter);
            } else if (ST_Tabular) {
                printTabular(*iter);
            }
            iter++;
        }
        if (ST_AggregateStats) {
            calculateAgregateSTats(&agregate_stats);
            printAggregate(&agregate_stats);
        }
    } catch (xercesc::DOMException& e ) {
        char * c_msg = tc(e.getMessage());
        std::cerr<<c_msg<<std::endl;
        xr(&c_msg);
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
    return 0;
}
void StatTool::parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
{
    
    StatManager * sm = new StatManager();
    ST_StatsVec.push_back(sm);
    char * c_cons = tc(parentNode->getAttribute(xmlParser.getDrseq()));
    std::string concensusRepeat =c_cons;
    sm->setConcensus(concensusRepeat);
    xr(&c_cons);

    char * c_gid = tc(parentNode->getAttribute(xmlParser.getGid()));
    std::string gid = c_gid;
    xr(&c_gid);
    sm->setGid(gid);
    
    
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getData())) {
            parseData(currentElement, xmlParser, sm);
        } /*else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getAssembly())) {
            if (ST_AssemblyStats) {
                parseAssembly(currentElement, xmlParser);
            }
        }*/
    }
}

void StatTool::parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getDrs())) {
            // change the direct repeats
            parseDrs(currentElement, xmlParser, statManager);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getSpacers())) {
            // change the spacers
            parseSpacers(currentElement, xmlParser, statManager);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getFlankers())) {
            // change the flankers
            parseFlankers(currentElement, xmlParser, statManager);
        }
    }
}

void StatTool::parseDrs(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {

        char * c_repeat = tc(currentElement->getAttribute(xmlParser.getSeq()));

        std::string repeat = c_repeat;

        xr(&c_repeat);
        statManager->addRepLenVec(static_cast<int>(repeat.length()));
        statManager->incrementRpeatCount();
    }
}

void StatTool::parseSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {

        char * c_spacer = tc(currentElement->getAttribute(xmlParser.getSeq()));
        std::string spacer = c_spacer;

        xr(&c_spacer);
        statManager->addSpLenVec(static_cast<int>(spacer.length()));
        char * c_cov = tc(currentElement->getAttribute(xmlParser.getCov()));
        std::string cov = c_cov;
        xr(&c_cov);
        if (!cov.empty()) {
            int cov_int;
            from_string(cov_int, cov, std::dec);
            statManager->addSpCovVec(cov_int);
            statManager->incrementSpacerCount();
        }
    }
}

void StatTool::parseFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager)
{
    
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        
        char * c_flanker = tc(currentElement->getAttribute(xmlParser.getSeq()));
        std::string flanker = c_flanker;
        xr(&c_flanker);
        statManager->addFlLenVec(static_cast<int>(flanker.length()));
        statManager->incrementFlankerCount();
    }
}

//void StatTool::parseAssembly(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
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
//            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getContig())) {
//                currentElement->setAttribute(xmlParser.getCid(), xmlParser.STR_2_XMLCH(getNextContigS()));
//                incrementContig();
//                parseContig(currentElement, xmlParser);
//            }
//            
//        }
//    }
//    
//}
//
//void StatTool::parseContig(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
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
//            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getCspacer())) {
//                if (ST_Spacers) {
//                    std::string spid = xmlParser.XMLCH_2_STR( currentElement->getAttribute(xmlParser.getSpid()));
//                    currentElement->setAttribute(xmlParser.getSpid(), xmlParser.STR_2_XMLCH(ST_SpacerMap[spid]));
//                }
//                if (ST_Spacers || ST_Repeats || ST_Flank) {
//                    parseCSpacer(currentElement, xmlParser);
//                }
//            }
//        }
//    }
//    
//}
//
//void StatTool::parseCSpacer(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
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
//            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getBspacers())) {
//                if (ST_Spacers || ST_Repeats) {
//                    parseLinkSpacers(currentElement, xmlParser);
//                }
//            } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getFspacers())) {
//                if (ST_Spacers || ST_Repeats) {
//                    parseLinkSpacers(currentElement, xmlParser);
//                }
//            } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getBflankers())) {
//                if (ST_Flank || ST_Repeats) {
//                    parseLinkFlankers(currentElement, xmlParser);
//                }
//            } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getFflankers())) {
//                if (ST_Flank || ST_Repeats) {
//                    parseLinkFlankers(currentElement, xmlParser);
//                }
//            }
//        }
//    }
//}
//
//void StatTool::parseLinkSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
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
//            if (ST_Spacers) {
//                std::string spid = xmlParser.XMLCH_2_STR( currentElement->getAttribute(xmlParser.getSpid()));
//                currentElement->setAttribute(xmlParser.getSpid(), xmlParser.STR_2_XMLCH( ST_SpacerMap[spid]));
//            }
//            if (ST_Repeats) {
//                std::string drid = xmlParser.XMLCH_2_STR( currentElement->getAttribute(xmlParser.getDrid()));
//                currentElement->setAttribute(xmlParser.getDrid(), xmlParser.STR_2_XMLCH( ST_RepeatMap[drid]));
//            }
//        }
//    }
//}
//
//void StatTool::parseLinkFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
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
//            if (ST_Flank) {
//                std::string flid = xmlParser.XMLCH_2_STR( currentElement->getAttribute(xmlParser.getFlid()));
//                currentElement->setAttribute(xmlParser.getFlid(), xmlParser.STR_2_XMLCH( ST_FlankMap[flid]));
//            }
//            
//            if (ST_Repeats) {
//                std::string drid = xmlParser.XMLCH_2_STR( currentElement->getAttribute(xmlParser.getDrid()));
//                currentElement->setAttribute(xmlParser.getDrid(), xmlParser.STR_2_XMLCH( ST_RepeatMap[drid]));
//            }
//        }
//    }
//}
void StatTool::calculateAgregateSTats(AStats * agregateStats)
{
    std::vector<StatManager * >::iterator iter;
    for(iter = begin(); iter != end(); iter++) {
        agregateStats->total_groups++;
        agregateStats->total_dr += (*iter)->getRpeatCount();
        agregateStats->total_dr_length += (*iter)->meanRepeatL();
        agregateStats->total_spacers += (*iter)->getSpacerCount();
        agregateStats->total_spacer_length += ((*iter)->getSpLenVec().empty()) ?  0 : (*iter)->meanSpacerL();
        agregateStats->total_spacer_cov +=((*iter)->getSpCovVec().empty()) ?  0 : (*iter)->meanSpacerC();
        agregateStats->total_flanker += (*iter)->getFlankerCount();
        agregateStats->total_flanker_length += ((*iter)->getFlLenVec().empty()) ? 0 : (*iter)->meanFlankerL();
    }
}
void StatTool::prettyPrint(StatManager * sm)
{
    std::cout<<sm->getGid()<<" | "<<sm->getConcensus()<<" | ";
    int i = 0;
    for ( i = 0; i < sm->getRpeatCount(); ++i) {
        std::cout<<REPEAT_CHAR;
    }
    for ( i = 0; i < sm->getSpacerCount() ; ++i) {
        std::cout<<SPACER_CHAR;
    }
    for ( i = 0; i < sm->getFlankerCount(); ++i) {
        std::cout<<FLANKER_CHAR;
    }
    std::cout<<"{ "<<sm->getRpeatCount()<< " " <<sm->getSpacerCount()<<" "<<sm->getFlankerCount()<<" } "<<std::endl;
}

void StatTool::printHeader()
{
    std::cout<<"GID"<<ST_Separator;
    std::cout<<"DR concensus"<<ST_Separator;
    std::cout<<"# DR Variants"<<ST_Separator;
    std::cout<<"Ave. DR Length"<<ST_Separator;
    std::cout<<"# spacers"<<ST_Separator;
    std::cout<<"Ave. SP Length"<<ST_Separator;
    std::cout<<"Ave. SP Cov"<<ST_Separator;
    std::cout<<"# Flankers"<<ST_Separator;
    std::cout<<"Ave. FL Length"<<std::endl;
    ST_WithHeader = false;
}

void StatTool::printTabular(StatManager * sm)
{
    if (ST_WithHeader) {
        printHeader();
    }
    std::cout<< sm->getGid()<<ST_Separator;
    std::cout<< sm->getConcensus()<<ST_Separator;
    std::cout<< sm->getRpeatCount()<< ST_Separator;
    std::cout<< sm->meanRepeatL()<<ST_Separator;
    std::cout<< sm->getSpacerCount()<<ST_Separator;
    if (!sm->getSpLenVec().empty()) {
        std::cout<< sm->meanSpacerL()<<ST_Separator;
    } else {
        std::cout<<0<<ST_Separator;
    }
    if (sm->getSpCovVec().empty()) {
        std::cout<<0<<ST_Separator;
    } else {
        std::cout<< sm->meanSpacerC()<<ST_Separator;
    }
    std::cout<< sm->getFlankerCount()<<ST_Separator;
    if (sm->getFlLenVec().empty()) {
        std::cout<<0<<std::endl;
    } else {
        std::cout<< sm->meanFlankerL()<<std::endl;
    }
}
void StatTool::printAggregate( AStats * agregate_stats)
{
    // if we pass through a NULL pointer then we print
    if (ST_WithHeader) {
        // the way this is called we should only print
        // the header only if tabular isn't set
        printHeader();
    }
    std::cout<<agregate_stats->total_groups<<ST_Separator;
    std::cout<<"*"<<ST_Separator;
    std::cout<<agregate_stats->total_dr<<ST_Separator;
    std::cout<<agregate_stats->total_dr_length/agregate_stats->total_groups<<ST_Separator;
    std::cout<<agregate_stats->total_spacers<<ST_Separator;
    std::cout<<agregate_stats->total_spacer_length/agregate_stats->total_groups<<ST_Separator;
    std::cout<<agregate_stats->total_spacer_cov/agregate_stats->total_groups<<ST_Separator;
    std::cout<<agregate_stats->total_flanker<<ST_Separator;
    std::cout<<agregate_stats->total_flanker_length/agregate_stats->total_groups<<std::endl;
}
int statMain (int argc, char ** argv)
{
    try {
		StatTool st;
		int opt_index = st.processOptions (argc, argv);
		if (opt_index >= argc) {
			throw crispr::input_exception("No input file provided" );
            
        } else {
			// get cracking and process that file
			return st.processInputFile(argv[opt_index]);
		}
	} catch(crispr::input_exception& re) {
        std::cerr<<re.what()<<std::endl;
        statUsage();
        return 1;
    } catch(crispr::exception& ce ) {
		std::cerr<<ce.what()<<std::endl;
		return 1;
	}
    return 0;
}
void statUsage(void)
{
    std::cout<<CRISPRTOOLS_PACKAGE_NAME<<" stat [-aghpst] [--header] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
    std::cout<<"-a                  print out aggregate summary, can be combined with -t -p"<<std::endl;
    std::cout<<"-h					print this handy help message"<<std::endl;
    std::cout<<"-H                  print out column headers in tabular output"<<std::endl;
	std::cout<<"-g INT[,n]          a comma separated list of group IDs that you would like to see stats for."<<std::endl;
    std::cout<<"-p                  pretty print"<<std::endl;
    std::cout<<"-s                  separator string for tabular output [default: '\t']"<<std::endl;
    std::cout<<"-t                  tabular output"<<std::endl;
}
