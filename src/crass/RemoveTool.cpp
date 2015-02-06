/*
 *  RemoveTool.cpp is part of the crisprtools project
 *  
 *  Created by Connor Skennerton on 22/12/11.
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

#include <iostream>
#include <cstdio>
#include <getopt.h>
#include "RemoveTool.h"
#include "Exception.h"
#include "StlExt.h"
#include "parser.h"
#include "config.h"
#include "Utils.h"

int removeMain(int argc, char ** argv)
{

    try {
        
        std::set<std::string> groups;
        std::string output_file;
        bool remove_files = false;
        int opt_index = processRemoveOptions(argc, argv, groups, output_file, remove_files);
        if(argc <= opt_index) {
            throw crispr::input_exception("Please specify an input file");
        }
        
        crispr::xml::parser xml_obj;

        xercesc::DOMDocument * xml_doc = xml_obj.setFileParser(argv[opt_index]);
        
        xercesc::DOMElement * root_elem = xml_doc->getDocumentElement();
        
        if( !root_elem ) throw(crispr::xml_exception(__FILE__, 
                                                     __LINE__, 
                                                     __PRETTY_FUNCTION__, 
                                                     "empty XML document" ));

        // get the children
        std::vector<xercesc::DOMElement * > bad_children;
        for (xercesc::DOMElement * currentElement = root_elem->getFirstElementChild(); 
             currentElement != NULL; 
             currentElement = currentElement->getNextElementSibling()) {
            if (xercesc::XMLString::equals(currentElement->getTagName(), xml_obj.tag_Group())) {
                // new group
                char * c_group_id = tc(currentElement->getAttribute(xml_obj.attr_Gid()));
                std::string group_id = c_group_id;
                if (groups.find(group_id.substr(1)) != groups.end() ) {
                    bad_children.push_back(currentElement);
                    if (remove_files) {
                        removeAssociatedData(currentElement, xml_obj);
                    }
                }
                xr(&c_group_id);
            }
        }
        std::vector<xercesc::DOMElement * >::iterator iter = bad_children.begin();
        while (iter != bad_children.end()) {
            root_elem->removeChild(*iter);
            iter++;
        }
        
        if (output_file.empty()) {
            xml_obj.printDOMToFile(argv[opt_index],xml_doc);
        } else {
            xml_obj.printDOMToFile(output_file);
        }
        
    } catch( xercesc::XMLException& e ) {
        char * message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        throw (crispr::xml_exception(__FILE__, 
                                     __LINE__,
                                     __PRETTY_FUNCTION__,
                                     (errBuf.str()).c_str()));
        xercesc::XMLString::release( &message );
    } catch (xercesc::DOMException& e) {
        char * message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        throw (crispr::xml_exception(__FILE__, 
                                     __LINE__,
                                     __PRETTY_FUNCTION__,
                                     (errBuf.str()).c_str()));
        xercesc::XMLString::release( &message );
    } catch (crispr::xml_exception& xe) {
        std::cerr<< xe.what()<<std::endl;
        return 1;
    } catch (crispr::input_exception& ie) {
        std::cerr<<ie.what()<<std::endl;
        removeUsage();
        return 1;
    }
    return 0;
}

void removeAssociatedData(xercesc::DOMElement * groupElement, 
                          crispr::xml::writer& xmlParser)
{
    for (xercesc::DOMElement * currentElement = groupElement->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_Metadata())) {
            parseMetadata(currentElement, xmlParser);
        }
    }
}

void parseMetadata(xercesc::DOMElement * parentNode, 
                   crispr::xml::writer& xmlParser) {
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); 
         currentElement != NULL; 
         currentElement = currentElement->getNextElementSibling()) {
        if (xercesc::XMLString::equals(currentElement->getTagName(), xmlParser.tag_File())) {
            char * c_url = tc(currentElement->getAttribute(xmlParser.attr_Url()));
            if (remove(c_url)) {
                perror("Cannot remove file");
            }
            xr(&c_url);
        }
    }
}

void removeUsage(void)
{
    std::cout<<PACKAGE_NAME<<" rm [-hor] -g <groups> file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
    std::cout<<"-h					print this handy help message"<<std::endl;
	std::cout<<"-g INT[,n]          a comma separated list of group IDs that you would like to remove"<<std::endl;
    std::cout<<"-o FILE             output file name. Default behaviour changes file inplace"<<std::endl;
    std::cout<<"-r                  Remove associated files"<<std::endl;
}
int processRemoveOptions(int argc, char ** argv, std::set<std::string>& groups, std::string& outputFile, bool& rem )
{
    try {
        int c;
        int index;
        static struct option long_options [] = {       
            {"help", no_argument, NULL, 'h'},
            {"groups",required_argument, NULL, 'g'},
            {"remove-file", no_argument, NULL, 'r'},
            {"outfile",required_argument,NULL, 'o'},
            {0,0,0,0}
        };
        while((c = getopt_long(argc, argv, "hg:o:r", long_options, &index)) != -1)
        {
            switch(c)
            {
                case 'h':
                {
                    removeUsage ();
                    exit(1);
                    break;
                }
                case 'g':
                {
                    if(fileOrString(optarg)) {
                        // its a file
                        parseFileForGroups(groups, optarg);
                        //ST_Subset = true;
                        
                    } else {
                        // its a string 
                        generateGroupsFromString(optarg, groups);
                    }
                    break;
                }
                    
                case 'o':
                {
                    outputFile = optarg;
                    break;
                } 
                case 'r':
                {
                    rem = true;
                }
                default:
                {
                    removeUsage();
                    exit(1);
                    break;
                }
            }
        }
        if (groups.empty()) 
        {
            throw crispr::input_exception("Please specify the groups to remove with -g");
        }
    } catch (crispr::input_exception& e) {
        std::cerr<<e.what()<<std::endl;
        removeUsage();
        exit(1);
    } catch (crispr::runtime_exception& e) {
        std::cerr<<e.what()<<std::endl;
        exit(1);
    }
    return optind;

}
