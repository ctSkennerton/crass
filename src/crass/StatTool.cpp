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
#include "config.h"
#include "Exception.h"
#include "reader.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
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

//void StatTool::generateGroupsFromString ( std::string str)
//{
//	std::set<std::string> vec;
//	split ( str, vec, ",");
//	ST_Groups = vec;
//	ST_Subset = true;
//}
int StatTool::processOptions (int argc, char ** argv)
{
	int c, index;
    struct option long_opts [] = { 
        {"header", no_argument, NULL, 'H'},
        {"coverage", no_argument, NULL, 0}
    };
	while((c = getopt_long(argc, argv, "ahHg:pPs:o:", long_opts, &index)) != -1)
	{
        switch(c)
		{
			case 'a':
            {
                ST_AggregateStats = true;
                break;
            }
            case 'p':
            {
                ST_OutputStyle = pretty;
                break;
            }
            case 'P':
            {
                ST_OutputStyle = veryPretty;
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

            case 'g':
            {
                if(fileOrString(optarg)) {
                    // its a file
                    parseFileForGroups(ST_Groups, optarg);
                    //ST_Subset = true;

                } else {
                    // its a string 
                    generateGroupsFromString(optarg, ST_Groups);
                }
                ST_Subset = true;
                break;
            }
            case 's':
            {
                ST_Separator = optarg;
                break;
            }

            case 'H':
            {
                ST_WithHeader = true;
                break;
            }
            case 0:
            {
                if (! strcmp("coverage", long_opts[index].name)) {
                    ST_DetailedCoverage = true;
                    ST_OutputStyle = coverage;
                }
                break;
            }
            default:
            {
                statUsage();
                exit(1);
                break;
            }
		}
	}
	return optind;
}

int StatTool::processInputFile(const char * inputFile)
{
    try {
        crispr::xml::reader xml_parser;
        std::ifstream in_file_stream(inputFile);
        if (in_file_stream.good()) {
            in_file_stream.close();
        } else {
            throw crispr::input_exception("cannot open input file");
        }
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        if (!root_elem) {
            throw crispr::xml_exception(__FILE__, 
                                        __LINE__, 
                                        __PRETTY_FUNCTION__, 
                                        "problem when parsing xml file");
        }
        int num_groups_to_process = static_cast<int>(ST_Groups.size());
        //std::cout<<num_groups_to_process<<std::endl;
        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild();
             currentElement != NULL; 
             currentElement = currentElement->getNextElementSibling()) {

            if (ST_Subset && num_groups_to_process == 0) {
                break;
            }
            // is this a group element
            if (xercesc::XMLString::equals(currentElement->getTagName(), xml_parser.tag_Group())) {
                char * c_gid = tc(currentElement->getAttribute(xml_parser.attr_Gid()));
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
        agregate_stats.total_reads = 0;
        // go through each of the groups and print out a pretty picture
        std::vector<StatManager *>::iterator iter = this->begin();
        int longest_consensus = 0;
        int longest_gid = 0;
        if (ST_OutputStyle == veryPretty) {
            while (iter != this->end()) {
                if(static_cast<int>((*iter)->getConcensus().length()) > longest_consensus) {
                    longest_consensus = static_cast<int>((*iter)->getConcensus().length());
                }
                if (static_cast<int>((*iter)->getGid().length()) > longest_gid) {
                    longest_gid = static_cast<int>((*iter)->getGid().length());
                }
            }
            iter = this->begin();
        }
        
        while (iter != this->end()) {
            switch (ST_OutputStyle) {
                case tabular:
                    printTabular(*iter);
                    break;
                case pretty:
                    prettyPrint(*iter);
                    break;
                case veryPretty:
                    veryPrettyPrint(*iter, longest_consensus, longest_gid);
                    break;
                case coverage:
                    printCoverage(*iter);
                    break;
                default:
                    break;
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
        return 1;
    }  catch (crispr::exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
    return 0;
}
void StatTool::parseGroup(xercesc::DOMElement * parentNode, 
                          crispr::xml::base& xmlParser)
{
    
    StatManager * sm = new StatManager();
    ST_StatsVec.push_back(sm);
    char * c_cons = tc(parentNode->getAttribute(xmlParser.attr_Drseq()));
    std::string concensusRepeat =c_cons;
    sm->setConcensus(concensusRepeat);
    xr(&c_cons);

    char * c_gid = tc(parentNode->getAttribute(xmlParser.attr_Gid()));
    std::string gid = c_gid;
    xr(&c_gid);
    sm->setGid(gid);
    
    
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Data())) {
            parseData(currentElement, xmlParser, sm);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Metadata())) {
            parseMetadata(currentElement, xmlParser, sm);
        }
    }
}

void StatTool::parseData(xercesc::DOMElement * parentNode, 
                         crispr::xml::base& xmlParser, 
                         StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Drs())) {
            // change the direct repeats
            parseDrs(currentElement, xmlParser, statManager);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacers())) {
            // change the spacers
            parseSpacers(currentElement, xmlParser, statManager);
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flankers())) {
            // change the flankers
            parseFlankers(currentElement, xmlParser, statManager);
        }
    }
}

void StatTool::parseDrs(xercesc::DOMElement * parentNode, 
                        crispr::xml::base& xmlParser, 
                        StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        char * c_repeat = tc(currentElement->getAttribute(xmlParser.attr_Seq()));

        std::string repeat = c_repeat;

        xr(&c_repeat);
        statManager->addRepLenVec(static_cast<int>(repeat.length()));
        statManager->incrementRpeatCount();
    }
}

void StatTool::parseSpacers(xercesc::DOMElement * parentNode, 
                            crispr::xml::base& xmlParser, 
                            StatManager * statManager)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        char * c_spacer = tc(currentElement->getAttribute(xmlParser.attr_Seq()));
        std::string spacer = c_spacer;
        xr(&c_spacer);
        statManager->addSpLenVec(static_cast<int>(spacer.length()));
        char * c_cov = tc(currentElement->getAttribute(xmlParser.attr_Cov()));
        std::string cov = c_cov;
        xr(&c_cov);
        if (!cov.empty()) {
            int cov_int;
            from_string(cov_int, cov, std::dec);
            statManager->addSpCovVec(cov_int);
        }
        statManager->incrementSpacerCount();
    }
}

void StatTool::parseFlankers(xercesc::DOMElement * parentNode, 
                             crispr::xml::base& xmlParser, 
                             StatManager * statManager)
{
    
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        char * c_flanker = tc(currentElement->getAttribute(xmlParser.attr_Seq()));
        std::string flanker = c_flanker;
        xr(&c_flanker);
        statManager->addFlLenVec(static_cast<int>(flanker.length()));
        statManager->incrementFlankerCount();
    }
}

void StatTool::parseMetadata(xercesc::DOMElement * parentNode, 
                             crispr::xml::base& xmlParser, 
                             StatManager * statManager) {
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_File())) {
            char * c_type_attr = tc(currentElement->getAttribute(xmlParser.attr_Type()));
            if (! strcmp(c_type_attr, "sequence")) {
                char * c_url = tc(currentElement->getAttribute(xmlParser.attr_Url()));
                statManager->setReadCount(calculateReads(c_url));
                xr(&c_url);
            }
            xr(&c_type_attr);
        }
    }
}
int StatTool::calculateReads(const char * fileName) {
    std::fstream sequence_file;
    sequence_file.open(fileName);
    int sequence_counter = 0;
    if (! sequence_file.good()) {
        
    }
    std::string line;
    while (sequence_file >> line) {
        if (line.substr(0,1) == ">") {
            sequence_counter++;
        }
    }
    return sequence_counter;
}

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
        agregateStats->total_reads += (*iter)->getReadCount();
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

void StatTool::veryPrettyPrint(StatManager * sm, int longestConsensus, int longestGID )
{
    //const int num_columns = 120;
    
    //int num_lines = sm->getSpacerCount() / num_columns;
    int current_gid_length = static_cast<int>(sm->getGid().length());
    int current_consensus_length = static_cast<int>(sm->getConcensus().length());
    int consensus_padding = longestConsensus - current_consensus_length;
    int gid_padding = longestGID - current_gid_length;
    
    //int total_length_before_spacers = current_gid_length + current_consensus_length + gid_padding + consensus_padding;
    
    std::cout<<sm->getGid();
    int i = 0;
    for ( i = 0; i < gid_padding; i++) {
        std::cout<<' ';
    }
    
    std::cout<<" | "<<sm->getConcensus();
    for ( i = 0; i < consensus_padding; i++) {
        std::cout<<' ';
    }
    std::cout<<" | ";
    
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
    std::cout<<"Ave. FL Length"<<ST_Separator;
    std::cout<<"# Reads"<<std::endl;
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
        std::cout<<0<<ST_Separator;
    } else {
        std::cout<< sm->meanFlankerL()<<ST_Separator;
    }
    std::cout<<sm->getReadCount()<<std::endl;
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
    if (agregate_stats->total_groups != 0) {
        std::cout<<agregate_stats->total_dr_length/agregate_stats->total_groups<<ST_Separator;
    } else {
        std::cout<<0<<ST_Separator;
    }
    std::cout<<agregate_stats->total_spacers<<ST_Separator;
    if (agregate_stats->total_groups != 0) {
        std::cout<<agregate_stats->total_spacer_length/agregate_stats->total_groups<<ST_Separator;
    } else {
        std::cout<<0<<ST_Separator;
    }
    if (agregate_stats->total_groups != 0) {
        std::cout<<agregate_stats->total_spacer_cov/agregate_stats->total_groups<<ST_Separator;
    } else {
        std::cout<<0<<ST_Separator;
    }
    std::cout<<agregate_stats->total_flanker<<ST_Separator;
    if (agregate_stats->total_groups != 0) {
        std::cout<<agregate_stats->total_flanker_length/agregate_stats->total_groups<<ST_Separator;
    } else {
        std::cout<<0<<ST_Separator;
    }
    if (agregate_stats->total_groups != 0) {
        std::cout<<agregate_stats->total_reads/agregate_stats->total_groups<<std::endl;
    } else {
        std::cout<<0<<std::endl;
    }
}

void StatTool::printCoverage(StatManager * sm)
{
    std::cout<< sm->getGid()<<ST_Separator;
    std::cout<< sm->getConcensus()<<ST_Separator;
    std::map<int, int> histogram;
    std::vector<int> coverage = sm->getSpCovVec();
    std::vector<int>::iterator iter;
    for (iter = coverage.begin(); iter != coverage.end(); iter++) {
        addOrIncrement(histogram, *iter);
    }
    std::map<int, int>::iterator h_iter;
    for (h_iter = histogram.begin(); h_iter != histogram.end(); h_iter++) {
        std::cout<<h_iter->first<<":"<<h_iter->second<<",";
    }
    std::cout<<std::endl;

    
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
    std::cout<<PACKAGE_NAME<<" stat [-aghpst] [--header] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
    std::cout<<"-a                  print out aggregate summary, can be combined with -t -p"<<std::endl;
    std::cout<<"-h					print this handy help message"<<std::endl;
    std::cout<<"-H                  print out column headers in tabular output"<<std::endl;
	std::cout<<"-g INT[,n]          a comma separated list of group IDs that you would like to see stats for."<<std::endl;
    std::cout<<"-p                  pretty print"<<std::endl;
    std::cout<<"-s                  separator string for tabular output [default: '\t']"<<std::endl;
    std::cout<<"-t                  tabular output"<<std::endl;
    std::cout<<"--coverage          Create a detailed report on the spacer coverage for each group"<<std::endl;
}
