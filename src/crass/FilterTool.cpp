// FilterTool.cpp
//
// Copyright (C) 2011, 2012 - Connor Skennerton
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

#include "FilterTool.h"
#include <libcrispr/Exception.h>
#include "config.h"
#include <libcrispr/writer.h>
#include <libcrispr/parser.h>
#include <libcrispr/StlExt.h>
#include <iostream>
#include <getopt.h>
#include "Utils.h"

#define exists(container,searchThing )  container.find(searchThing) != container.end()


int FilterTool::processOptions (int argc, char ** argv)
{
	int c;
    int index;
    static struct option long_options [] = {       
        {"help", no_argument, NULL, 'h'},
        {"outfile", required_argument, NULL, 'o'},
        {"spacer",required_argument,NULL,'s'},
        {"direct-repeat", required_argument, NULL, 'd'},
        {"flanker", required_argument, NULL, 'f'},
        {"coverage",required_argument,NULL,'C'},
        {0,0,0,0}
    };
	while((c = getopt_long(argc, argv, "hs:c:f:d:o:C:", long_options,&index)) != -1)
	{
        switch(c)
		{
            case 'h':
			{
				filterUsage();
				exit(1);
				break;
			}
			case 's':
			{
				if (! from_string(FT_Spacers, optarg, std::dec)) {
                    crispr::input_exception("cannot convert arguement of -s to integer");
                }
                break;
			}
            case 'o':
            {
                FT_OutputFile = optarg;
                break;
            }
            case 'f':
            {
                if (! from_string(FT_Flank, optarg, std::dec)) {
                    crispr::input_exception("cannot convert arguement of -f to integer");
                }
                break;
            }
            case 'd':
            {
                if(! from_string(FT_Repeats, optarg, std::dec)){
                    crispr::input_exception("cannot convert arguement of -d to integer");
                }
                break;
            }
            case 'c':
            {
                if(! from_string(FT_contigs, optarg, std::dec)){
                    crispr::input_exception("cannot convert arguement of -c to integer");
                }
                break;
            }
            case 'C':
            {
                if(! from_string(FT_Coverage, optarg, std::dec)){
                    crispr::input_exception("cannot convert arguement of -C to integer");
                }
                break;
            }
            default:
            {
                filterUsage();
                exit(1);
                break;
            }
		}
	}
	return optind;
}

int FilterTool::processInputFile(const char * inputFile)
{
    try {
        crispr::xml::writer output_xml;
        int output_error;
        xercesc::DOMElement * output_root_elem = output_xml.createDOMDocument("crispr", "1.1", output_error);
        if(output_error && NULL == output_root_elem) {
            throw crispr::xml_exception(__FILE__, 
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        "Cannot create output xml file");
        }
        
        crispr::xml::parser xml_parser;
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        if (!root_elem) {
            throw crispr::xml_exception(__FILE__, 
                                        __LINE__, 
                                        __PRETTY_FUNCTION__, 
                                        "problem when parsing xml file");
        }
        if (FT_OutputFile.empty()) {
            FT_OutputFile = inputFile;
        }
        
        std::vector<xercesc::DOMElement * > good_children;

        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); 
             currentElement != NULL; 
             currentElement = currentElement->getNextElementSibling()) {
                
            // the user wants to change any of these 
            if (FT_Spacers || FT_Repeats || FT_Flank || FT_Coverage) {
                if (! parseGroup(currentElement, xml_parser)) {
                    good_children.push_back(currentElement);
                }
            }
        }
        std::vector<xercesc::DOMElement * >::iterator iter = good_children.begin();
        xercesc::DOMDocument * output_document = output_xml.getDocumentObj();
        while (iter != good_children.end()) {
            output_root_elem->appendChild(output_document->importNode(*iter, true));
            iter++;
        }

        output_xml.printDOMToFile(FT_OutputFile, output_document);
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    } catch (xercesc::DOMException& e) {
        char * c_msg = tc(e.getMessage());
        XERCES_STD_QUALIFIER cerr << "An error occurred during creation of output transcoder. Msg is:"
        << XERCES_STD_QUALIFIER endl
        << c_msg << XERCES_STD_QUALIFIER endl;
        xr(&c_msg);
        return 1;
    }
    
    return 0;
}

// return true if group should be removed
bool FilterTool::parseGroup(xercesc::DOMElement * parentNode, 
                            crispr::xml::reader& xmlParser)
{
    // get the data tag and make sure that everything is good
    xercesc::DOMElement * currentElement = parentNode->getFirstElementChild();
    if (NULL != currentElement) {
        std::set<std::string> spacers_to_remove;

        if( parseData(currentElement, xmlParser, spacers_to_remove) ) return true;  
        
        // get the assembly node
        currentElement = parentNode->getLastElementChild();
        parseAssembly(currentElement, xmlParser, spacers_to_remove);
    } else {
        throw crispr::xml_exception(__FILE__,
                                    __LINE__,
                                    __PRETTY_FUNCTION__,
                                    "There is no Data");
    }
    return false;
}


// return true if group should be removed
bool FilterTool::parseData(xercesc::DOMElement * parentNode, 
                           crispr::xml::reader& xmlParser,
                           std::set<std::string>& spacersToRemove)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Drs())) {
            if (FT_Repeats) {
                // change the direct repeats
                if (FT_Repeats > parseDrs(currentElement)) {
                    return true;
                }
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacers())) {
            if (FT_Spacers | FT_Coverage) {
                // change the spacers
                if (FT_Spacers > parseSpacers(currentElement, xmlParser, spacersToRemove)) {
                    return true;
                }
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flankers())) {
            if (FT_Flank) {
                // change the flankers
                if ( FT_Flank > parseFlankers(currentElement)) {
                    return true;
                }
            }
        }
    }
    return false;
}
int FilterTool::countElements(xercesc::DOMElement * parentNode)
{
    int count = 0; 
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        ++count;
    }
    return count;
}

int FilterTool::parseSpacers(xercesc::DOMElement *parentNode, crispr::xml::reader& xmlParser, std::set<std::string>& spacersToRemove) {
    if (FT_Coverage) {
        std::vector<xercesc::DOMElement * > remove_list;
        for (xercesc::DOMElement * currentSpacer = parentNode->getFirstElementChild(); 
             currentSpacer != NULL; 
             currentSpacer = currentSpacer->getNextElementSibling()) {
            char * c_spacer_coverage = tc(currentSpacer->getAttribute(xmlParser.attr_Cov()));
            int cov;
            from_string(cov, c_spacer_coverage, std::dec);
            xr(&c_spacer_coverage);
            if (cov < FT_Coverage) {
                // remove spacer
                remove_list.push_back(currentSpacer);
                char * c_spid = tc(currentSpacer->getAttribute(xmlParser.attr_Spid()));
                spacersToRemove.insert(c_spid);
                xr(&c_spid);
            }
        }
        std::vector<xercesc::DOMElement * >::iterator iter;
        for (iter = remove_list.begin(); iter != remove_list.end(); iter++) {
            parentNode->removeChild(*iter);
        }
    } else if(FT_Spacers) {
        return countElements(parentNode);
    }
    return 1;
}

void FilterTool::parseAssembly(xercesc::DOMElement * parentNode, 
                             crispr::xml::reader& xmlParser, std::set<std::string>& spacersToRemove)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Contig())) {
            char * c_contig_id = tc(currentElement->getAttribute(xmlParser.attr_Cid()));
            std::string contig_id = c_contig_id;
            xr(&c_contig_id);
            parseContig(currentElement, xmlParser, contig_id, spacersToRemove);
        }
        
        
    }
    
}

void FilterTool::parseContig(xercesc::DOMElement * parentNode, 
                           crispr::xml::reader& xmlParser, 
                             std::string& contigId, std::set<std::string>& spacersToRemove)
{
    std::vector<xercesc::DOMElement* > remove_list;
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Cspacer())) {
            // get the node
            char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            
            if(exists(spacersToRemove, c_spid)) {
                remove_list.push_back(currentElement);
                continue;
            }            
			xr(&c_spid);            
            parseCSpacer(currentElement, xmlParser, contigId, spacersToRemove);
        }
    }
    std::vector<xercesc::DOMElement* >::iterator iter;
    for (iter = remove_list.begin(); iter != remove_list.end(); iter++) {
        parentNode->removeChild(*iter);
    }
    
}

void FilterTool::parseCSpacer(xercesc::DOMElement * parentNode, 
                            crispr::xml::reader& xmlParser, 
                            std::string& contigId, std::set<std::string>& spacersToRemove)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Fspacers())) {
            parseLinkSpacers(currentElement, xmlParser, contigId, spacersToRemove);
            
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Bspacers())) {
            parseLinkSpacers(currentElement, xmlParser, contigId, spacersToRemove);
            
        }
        
    }
}

void FilterTool::parseLinkSpacers(xercesc::DOMElement * parentNode, 
                                crispr::xml::reader& xmlParser, 
                                std::string& contigId, std::set<std::string>& spacersToRemove)
{
    std::vector<xercesc::DOMElement* > remove_list;

    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        
        char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
        if(exists(spacersToRemove, c_spid)) {
            remove_list.push_back(currentElement);
            continue;
        } 
        xr(&c_spid);
    }
    std::vector<xercesc::DOMElement* >::iterator iter;
    for (iter = remove_list.begin(); iter != remove_list.end(); iter++) {
        parentNode->removeChild(*iter);
    }
}

int filterMain (int argc, char ** argv)
{
    try {
		FilterTool ft;
		int opt_index = ft.processOptions (argc, argv);
		if (opt_index >= argc) {
			throw crispr::input_exception("No input file provided" );
            
        } else {
			// get cracking and process that file
			return ft.processInputFile(argv[opt_index]);
		}
	} catch(crispr::input_exception& re) {
        std::cerr<<re.what()<<std::endl;
        filterUsage();
        return 1;
    } catch(crispr::exception& ce ) {
		std::cerr<<ce.what()<<std::endl;
		return 1;
	}
    return 0;
}

void filterUsage (void)
{
    std::cout<<PACKAGE_NAME<<" filter [-ohsdf] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h                  Print this handy help message"<<std::endl;
    std::cout<<"-o FILE             Output file name, creates a filtered copy of the input file  [default: modify input file inplace]" <<std::endl; 
	std::cout<<"-s INT              Filter based on the number of spacers the spacers "<<std::endl;
	std::cout<<"-d INT              Filter based on the direct repeats "<<std::endl;
	std::cout<<"-f INT              Filter based on the flanking sequences "<<std::endl;
    std::cout<<"-C INT              Filter based on spacer coverage"<<std::endl;
	//std::cout<<"-c INT				filter based on the contigs "<<std::endl;
}
