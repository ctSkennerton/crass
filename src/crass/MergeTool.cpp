// MergeTool.cpp
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

#include "MergeTool.h"
#include "Exception.h"
#include "writer.h"
#include "reader.h"
#include "config.h"
#include <getopt.h>
#include <sstream>

int MergeTool::processOptions (int argc, char ** argv)
{
	int c;
    int index;
    static struct option long_options [] = {       
        {"help", no_argument, NULL, 'h'},
        {"sanitise", no_argument, NULL, 's'},
        {"outfile", required_argument, NULL, 'o'},
        {0,0,0,0}

    };
	while((c = getopt_long(argc, argv, "hso:", long_options, &index)) != -1)
	{
        switch(c)
		{
			case 'h':
			{
				mergeUsage ();
				exit(0);
				break;
			}
			case 's':
			{
				MT_Sanitise = true;
                break;
			}
            case 'o':
            {
                MT_OutFile = optarg;
                break;
            }
            default:
            {
                mergeUsage();
                exit(1);
                break;
            }
		}
	}
	return optind;
}
int mergeMain (int argc, char ** argv)
{
	try {
        MergeTool mt;
        int opt_index = mt.processOptions(argc, argv);
        if (opt_index >= argc) {
			throw crispr::input_exception("No input files provided" );
            
        } else if (opt_index == argc - 1) {
			// less than 2 input files
			throw crispr::input_exception("You must provide at least two input files to merge");
		} else {
            // merge!!
            
            // create a DOM document
            crispr::xml::writer master_DOM;
            int master_DOM_error;
            xercesc::DOMElement * master_root_elem = master_DOM.createDOMDocument("crispr", "1.1", master_DOM_error);
            if (master_root_elem != NULL && master_DOM_error == 0) {
                
                xercesc::DOMDocument * master_doc = master_DOM.getDocumentObj();
                while (opt_index < argc) {
                    
                    // create a file parser
                    crispr::xml::reader input_file;                    
                    xercesc::DOMDocument * input_doc = input_file.setFileParser(argv[opt_index]);
                    // Get the top-level element: 
                    xercesc::DOMElement* elementRoot = input_doc->getDocumentElement();
                    
                    if( !elementRoot ) throw(crispr::xml_exception( __FILE__,
                                                                   __LINE__,
                                                                   __PRETTY_FUNCTION__,
                                                                   "empty XML document" ));
                    
                    // get the children
                    for (xercesc::DOMElement * currentElement = elementRoot->getFirstElementChild(); 
                         currentElement != NULL; 
                         currentElement = currentElement->getNextElementSibling()) {

                        if( xercesc::XMLString::equals(currentElement->getTagName(), input_file.tag_Group())) {
                            
                            if (mt.getSanitise()) {
                                // change the name 
                                std::stringstream ss;
                                ss <<'G'<< mt.getNextGroupID();
                                XMLCh * x_group = tc(ss.str().c_str());
                                currentElement->setAttribute(master_DOM.attr_Gid(), x_group);
                                mt.incrementGroupID();
                                xr(&x_group);
                            
                            } else {
                                // check if we already seen it if so warn the user
                                char * gid = tc(currentElement->getAttribute(master_DOM.attr_Gid()));
                                
                                if ( mt.find(gid) != mt.end()) {
                                    // this group id has been seen before
                                    std::cout<<"Group IDs in the two files conflict "<<gid<<" seen more than once."<<std::endl;
                                    std::cout<<"Try using -s to avoid this or use "<<PACKAGE_NAME<<" sanitise to fix these conflicts"<<std::endl;
                                    
                                } else {
                                    // add in to the set
                                    mt.insert(gid);
                                }
                            }
                            master_root_elem->appendChild(master_doc->importNode(currentElement, true));
                        }
                        
                    }
                    opt_index++;
                }   
                master_DOM.printDOMToFile(mt.getFileName());
            } else {
                throw crispr::xml_exception(__FILE__, 
                                            __LINE__,
                                            __PRETTY_FUNCTION__,
                                            "no root for master DOM");
            }

        }
    } catch (crispr::input_exception& e) {
        std::cerr<<e.what()<<std::endl;
        mergeUsage();
        return 1;
    } catch (crispr::xml_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 2;
    }
    
    return 0;
}

void mergeUsage(void)
{
	std::cout<<PACKAGE_NAME<<" merge [-hso] file1.crispr file2.crispr [1,n]"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h					print this handy help message"<<std::endl;
    std::cout<<"-o FILE             output file  [default: crisprtools_merged.crispr]" <<std::endl; 
	std::cout<<"-s					sanitise the names so that the resulting output file contains completely unique group IDs"<<std::endl;
}
