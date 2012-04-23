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
#include "Exception.h"
#include "config.h"
#include "StlExt.h"
#include <iostream>

int FilterTool::processOptions (int argc, char ** argv)
{
	int c;
	while((c = getopt(argc, argv, "hs:c:f:d:o:")) != -1)
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
        crispr::XML xml_parser;
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        if (!root_elem) {
            throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, "problem when parsing xml file");
        }
        if (FT_OutputFile.empty()) {
            FT_OutputFile = inputFile;
        }
        
        std::vector<xercesc::DOMElement * > bad_children;

        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
                
            // the user wants to change any of these 
            if (FT_Spacers || FT_Repeats || FT_Flank || FT_contigs) {
                if (parseGroup(currentElement, xml_parser)) {
                    bad_children.push_back(currentElement);
                }
            }
        }
        std::vector<xercesc::DOMElement * >::iterator iter = bad_children.begin();
        while (iter != bad_children.end()) {
            root_elem->removeChild(*iter);
            iter++;
        }

        xml_parser.printDOMToFile(FT_OutputFile, input_doc_obj);
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
    
    return 0;
}

// return true if group should be removed
bool FilterTool::parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
{
    // get the data tag and make sure that everything is good
    xercesc::DOMElement * currentElement = parentNode->getFirstElementChild();
    if (NULL != currentElement) {
        if (FT_Spacers || FT_Repeats || FT_Flank) {
            return parseData(currentElement, xmlParser);             
        }
    } 
    return false;
}


// return true if group should be removed
bool FilterTool::parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getDrs())) {
            if (FT_Repeats) {
                // change the direct repeats
                if (FT_Repeats > parseDrs(currentElement)) {
                    return true;
                }
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getSpacers())) {
            if (FT_Spacers) {
                // change the spacers
                if (FT_Spacers > parseSpacers(currentElement)) {
                    return true;
                }
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.getFlankers())) {
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
int FilterTool::parseDrs(xercesc::DOMElement * parentNode)
{
    int count = 0; 
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        
        ++count;
    }
    return count;
}
int FilterTool::parseSpacers(xercesc::DOMElement * parentNode)
{
    int count = 0;
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        
        ++count;
    }
    return count;
}

int FilterTool::parseFlankers(xercesc::DOMElement * parentNode)
{
    int count = 0;
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
        
        ++count;
        
    }
    return count;
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
    std::cout<<CRISPRTOOLS_PACKAGE_NAME<<" filter [-ohsdf] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h                  Print this handy help message"<<std::endl;
    std::cout<<"-o FILE             Output file name, creates a filtered copy of the input file  [default: modify input file inplace]" <<std::endl; 
	std::cout<<"-s INT              Filter based on the number of spacers the spacers "<<std::endl;
	std::cout<<"-d INT              Filter based on the direct repeats "<<std::endl;
	std::cout<<"-f INT              Filter based on the flanking sequences "<<std::endl;
	//std::cout<<"-c INT				filter based on the contigs "<<std::endl;
}
