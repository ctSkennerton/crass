// SanitiseTool.cpp
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

#include <getopt.h>
#include "SanitiseTool.h"
#include "Exception.h"
#include "parser.h"
#include "config.h"

#include <iostream>
int SanitiseTool::processOptions (int argc, char ** argv)
{
	int c;
    int index;
    static struct option long_options [] = {       
        {"help", no_argument, NULL, 'h'},
        {"all", no_argument, NULL, 'a'},
        //{"groups",required_argument, NULL, 'g'},
        {"spacer",no_argument,NULL,'s'},
        {"direct-repeat", no_argument, NULL, 'd'},
        {"flanker", no_argument, NULL, 'f'},
        {"contig", no_argument, NULL, 'c'},
        {"outfile",required_argument,NULL, 'o'},
        //{"outfile-dir",required_argument,NULL,'O'},
        {0,0,0,0}
    };
	while((c = getopt_long(argc, argv, "ahscfdo:", long_options, &index)) != -1)
	{
        switch(c)
		{
			case 'a':
            {
                ST_contigs = ST_Repeats = ST_Flank = ST_Spacers = true;
                break;
            }
            case 'h':
			{
				sanitiseUsage();
				exit(0);
				break;
			}
			case 's':
			{
				ST_Spacers = true;
                break;
			}
            case 'o':
            {
                ST_OutputFile = optarg;
                break;
            }
            case 'f':
            {
                ST_Flank = true;
                break;
            }
            case 'd':
            {
                ST_Repeats = true;
                break;
            }
            case 'c':
            {
                ST_contigs = true;
                break;
            }
            default:
            {
                sanitiseUsage();
                exit(1);
                break;
            }
		}
	}
    if (!(ST_contigs | ST_Spacers | ST_Repeats | ST_Flank)) {
        throw crispr::input_exception("Please specify one of -s -f -d -c");
    }
	return optind;
}

int SanitiseTool::processInputFile(const char * inputFile)
{
    try {
        crispr::xml::parser xml_parser;
        xercesc::DOMDocument * input_doc_obj = xml_parser.setFileParser(inputFile);
        xercesc::DOMElement * root_elem = input_doc_obj->getDocumentElement();
        if (!root_elem) {
            throw crispr::xml_exception(__FILE__, 
                                        __LINE__, 
                                        __PRETTY_FUNCTION__, 
                                        "problem when parsing xml file");
        }
        if (ST_OutputFile.empty()) {
            ST_OutputFile = inputFile;
        }

        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); 
             currentElement != NULL; 
             currentElement = currentElement->getNextElementSibling()) {

                
                // is this a group element
                if (xercesc::XMLString::equals(currentElement->getTagName(), xml_parser.tag_Group())) {
                    
                    XMLCh * x_next_group_num = tc( getNextGroupS().c_str());
                    currentElement->setAttribute(xml_parser.attr_Gid(), x_next_group_num);
                    incrementGroup();
                    xr(&x_next_group_num);
                    // the user wants to change any of these 
                    if (ST_Spacers || ST_Repeats || ST_Flank || ST_contigs) {
                        parseGroup(currentElement, xml_parser);
                    }
                }
            
            setNextRepeat(1);
            setNextContig(1);
            setNextRepeat(1);
            setNextSpacer(1);
        }
        xml_parser.printDOMToFile(ST_OutputFile, input_doc_obj);
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
    
    return 0;
}
void SanitiseTool::parseGroup(xercesc::DOMElement * parentNode, 
                              crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Data())) {
            if (ST_Spacers || ST_Repeats || ST_Flank) {
                parseData(currentElement, xmlParser);
            }
            
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Assembly())) {
            if (ST_contigs || ST_Spacers || ST_Repeats || ST_Flank) {
                parseAssembly(currentElement, xmlParser);
            }
        }
    }
}

void SanitiseTool::parseData(xercesc::DOMElement * parentNode, 
                             crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Drs())) {
            if (ST_Repeats) {
                // change the direct repeats
                parseDrs(currentElement, xmlParser);
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacers())) {
            if (ST_Spacers) {
                // change the spacers
                parseSpacers(currentElement, xmlParser);
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flankers())) {
            if (ST_Flank) {
                // change the flankers
                parseFlankers(currentElement, xmlParser);
            }
        }
    }
}
void SanitiseTool::parseDrs(xercesc::DOMElement * parentNode, 
                            crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Dr())) {
            char * c_drid = tc(currentElement->getAttribute(xmlParser.attr_Drid()));
            std::string drid = c_drid;
            ST_RepeatMap[drid] = getNextRepeatS();
            xr(&c_drid);
            XMLCh * x_next_repeat_num = tc(getNextRepeatS().c_str());
            currentElement->setAttribute(xmlParser.attr_Drid(), x_next_repeat_num);
            xr(&x_next_repeat_num);
            incrementRepeat();
        }
        
    }
}
void SanitiseTool::parseSpacers(xercesc::DOMElement * parentNode, 
                                crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Spacer())) {
            char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            std::string spid = c_spid;
            ST_SpacerMap[spid] = getNextSpacerS();
            xr(&c_spid);
            XMLCh * x_next_spacer_num = tc(getNextSpacerS().c_str());
            currentElement->setAttribute(xmlParser.attr_Spid(), x_next_spacer_num);
            xr(&x_next_spacer_num);
            incrementSpacer();
        }
    }
}

void SanitiseTool::parseFlankers(xercesc::DOMElement * parentNode, 
                                 crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Flanker())) {
            char * c_flid = tc(currentElement->getAttribute(xmlParser.attr_Flid()));
            std::string flid = c_flid;
            ST_FlankMap[flid] = getNextFlankerS();
            xr(&c_flid);
            XMLCh * x_next_flanker_num = tc(getNextFlankerS().c_str());
            currentElement->setAttribute(xmlParser.attr_Flid(), x_next_flanker_num);
            xr(&x_next_flanker_num);
            incrementFlanker();
        }
        
    }
}

void SanitiseTool::parseAssembly(xercesc::DOMElement * parentNode, 
                                 crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Contig())) {
            XMLCh * x_next_contig = tc(getNextContigS().c_str());
            currentElement->setAttribute(xmlParser.attr_Cid(), x_next_contig);
            incrementContig();
            xr(&x_next_contig);
            parseContig(currentElement, xmlParser);
        }
    }
}

void SanitiseTool::parseContig(xercesc::DOMElement * parentNode, 
                               crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Cspacer())) {
            if (ST_Spacers) {
                char * c_spid = tc( currentElement->getAttribute(xmlParser.attr_Spid()));
                std::string spid = c_spid;
                XMLCh * x_new_spid = tc(ST_SpacerMap[spid].c_str());
                currentElement->setAttribute(xmlParser.attr_Spid(), x_new_spid);
                xr(&c_spid);
                xr(&x_new_spid);
            }
            if (ST_Spacers || ST_Repeats || ST_Flank) {
                parseCSpacer(currentElement, xmlParser);
            }
        }
    }
}

void SanitiseTool::parseCSpacer(xercesc::DOMElement * parentNode, 
                                crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Bspacers())) {
            if (ST_Spacers || ST_Repeats) {
                parseLinkSpacers(currentElement, xmlParser);
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Fspacers())) {
            if (ST_Spacers || ST_Repeats) {
                parseLinkSpacers(currentElement, xmlParser);
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Bflankers())) {
            if (ST_Flank || ST_Repeats) {
                parseLinkFlankers(currentElement, xmlParser);
            }
        } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Fflankers())) {
            if (ST_Flank || ST_Repeats) {
                parseLinkFlankers(currentElement, xmlParser);
            }
        }
    }
}

void SanitiseTool::parseLinkSpacers(xercesc::DOMElement * parentNode, 
                                    crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (ST_Spacers) {
            char * c_spid = tc(currentElement->getAttribute(xmlParser.attr_Spid()));
            std::string spid = c_spid;
            XMLCh * x_new_spid = tc( ST_SpacerMap[spid].c_str());
            currentElement->setAttribute(xmlParser.attr_Spid(), x_new_spid);
            xr(&c_spid);
            xr(&x_new_spid);
        }
        if (ST_Repeats) {
            char * c_drid = tc( currentElement->getAttribute(xmlParser.attr_Drid()));
            std::string drid = c_drid;
            XMLCh * x_new_drid = tc(ST_RepeatMap[drid].c_str());
            currentElement->setAttribute(xmlParser.attr_Drid(), x_new_drid);
            xr(&c_drid);
            xr(&x_new_drid);
        }
        
    }
}

void SanitiseTool::parseLinkFlankers(xercesc::DOMElement * parentNode, 
                                     crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {

        if (ST_Flank) {
            char * c_flid = tc( currentElement->getAttribute(xmlParser.attr_Flid()));
            std::string flid = c_flid;
            XMLCh * x_new_flid = tc( ST_FlankMap[flid].c_str());
            currentElement->setAttribute(xmlParser.attr_Flid(), x_new_flid);
            xr(&c_flid);
            xr(&x_new_flid);
        }

        if (ST_Repeats) {
            char * c_drid = tc( currentElement->getAttribute(xmlParser.attr_Drid()));
            std::string drid = c_drid;
            XMLCh * x_new_drid = tc(ST_RepeatMap[drid].c_str());
            currentElement->setAttribute(xmlParser.attr_Drid(), x_new_drid);
            xr(&c_drid);
            xr(&x_new_drid);
        }
        
    }
}

int sanitiseMain (int argc, char ** argv)
{
    try {
		SanitiseTool st;
		int opt_index = st.processOptions (argc, argv);
		if (opt_index >= argc) {
			throw crispr::input_exception("No input file provided" );
            
        } else {
			// get cracking and process that file
			return st.processInputFile(argv[opt_index]);
		}
	} catch(crispr::input_exception& re) {
        std::cerr<<re.what()<<std::endl;
        sanitiseUsage();
        return 1;
    } catch(crispr::exception& ce ) {
		std::cerr<<ce.what()<<std::endl;
		return 1;
	}
}

void sanitiseUsage(void)
{
    std::cout<<PACKAGE_NAME<<" sanitise [-ohcsdfa] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h                  Print this handy help message"<<std::endl;
    std::cout<<"-o FILE             Output file name, creates a sanitised copy of the input file  [default: sanitise input file inplace]" <<std::endl; 
    std::cout<<"-a                  Sanitise all features. Equivelent to -sdfc"<<std::endl;
    std::cout<<"-s                  Sanitise the spacers "<<std::endl;
	std::cout<<"-d                  Sanitise the direct repeats "<<std::endl;
	std::cout<<"-f                  Sanitise the flanking sequences "<<std::endl;
    std::cout<<"-c                  Sanitise the contigs "<<std::endl;

    

}
