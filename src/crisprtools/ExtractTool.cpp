// ExtractTool.cpp
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
#include "ExtractTool.h"
#include "Utils.h"
#include "config.h"
#include "Exception.h"
#include "XML.h"

#include "StlExt.h"
#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>


ExtractTool::ExtractTool (void)
{
	ET_Spacer = false;
	ET_DirectRepeat = false;
	ET_Flanker = false;
	ET_SplitGroup = false;
	ET_SplitType = false;
	ET_Subset = false;
    ET_OutputPrefix = "./";
}

ExtractTool::~ExtractTool (void)
{}

void ExtractTool::generateGroupsFromString ( std::string str)
{
	std::set<std::string> vec;
	split ( str, vec, ",");
	ET_Group = vec;
	ET_Subset = true;
}

int ExtractTool::processOptions (int argc, char ** argv)
{
	int c;
    int index;
    static struct option long_options [] = {       
        {"help", required_argument, NULL, 'h'}
    };
	while((c = getopt_long(argc, argv, "hg:sdfxyo:", long_options, &index)) != -1)
	{
        switch(c)
		{
			case 'h':
			{
				extractUsage ();
				exit(1);
				break;
			}
			case 'g':
			{
				generateGroupsFromString (optarg);
				break;
			}
			case 's':
			{
				ET_Spacer = true;
                break;
			}
			case 'd':
			{
				ET_DirectRepeat = true;
                break;
			}
			case 'f':
			{
				ET_Flanker = true;
                break;
			}
			case 'x':
			{
				ET_SplitGroup = true;
                break;
			}
			case 'y':
			{
				ET_SplitType = true;
                break;
			}
            case 'o':
            {
                ET_OutputPrefix = optarg;
                // just in case the user put '.' or '..' or '~' as the output directory
                if (ET_OutputPrefix[ET_OutputPrefix.length() - 1] != '/')
                {
                    ET_OutputPrefix += '/';
                }
                
                // check if our output folder exists
                struct stat file_stats;
                if (0 != stat(ET_OutputPrefix.c_str(),&file_stats)) 
                {
                    recursiveMkdir(ET_OutputPrefix);
                }
                break;
            }
            default:
            {
                extractUsage();
                exit(1);
                break;
            }
		}
	}
	return optind;
}

int ExtractTool::processInputFile(const char * inputFile)
{
    // open the file
    crispr::XML xml_obj;
    try {
        xercesc::DOMDocument * xml_doc = xml_obj.setFileParser(inputFile);
        
        xercesc::DOMElement * root_elem = xml_doc->getDocumentElement();
        
        if( !root_elem ) throw(crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, "empty XML document" ));
        
        // neither -x or -y
        if ((!ET_SplitGroup) && (!ET_SplitType)) {
            // just open a single file handle
            ET_OneStream.open((ET_OutputPrefix +"extracted_data.fa").c_str());
        
        } else if(!ET_SplitGroup) {
            // -y but not -x
            // we are only spliting on type make three generic files
            if (ET_DirectRepeat) {
                ET_RepeatStream.open((ET_OutputPrefix + "direct_repeats.fa").c_str());
            }
            if (ET_Spacer) {
                ET_SpacerStream.open((ET_OutputPrefix + "spacers.fa").c_str());
            }
            if (ET_Flanker) {
                ET_FlankerStream.open((ET_OutputPrefix + "flankers.fa").c_str());
            }
        }
        
        parseWantedGroups(xml_obj, root_elem);
        
    } catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        throw (crispr::xml_exception(__FILE__, __LINE__,__PRETTY_FUNCTION__,(errBuf.str()).c_str()));
        xercesc::XMLString::release( &message );
        
    } catch (crispr::xml_exception& xe) {
        std::cerr<< xe.what()<<std::endl;
        return 1;
    }
    if ((!ET_SplitGroup) && (!ET_SplitType)) {
        // just open a single file handle
        ET_OneStream.close();
        
    } else if(!ET_SplitGroup) {
        // we are only spliting on type make three generic files
        if (ET_Spacer) {
            ET_SpacerStream.close();
        }
        if (ET_DirectRepeat) {
            ET_RepeatStream.close();
        }
        if (ET_Flanker) {
            ET_FlankerStream.close();
        }
    }

    return 0;
}

void ExtractTool::parseWantedGroups(crispr::XML& xmlObj, xercesc::DOMElement * rootElement)
{
    
    try {
        
        for (xercesc::DOMElement * currentElement = rootElement->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
            // is this a group element
            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlObj.getGroup())) {
                // new group
                char * c_group_id = tc(currentElement->getAttribute(xmlObj.getGid()));
                std::string group_id = c_group_id;
                if (ET_Subset) {
                    
                    // we only want some of the groups look at ET_Groups
                    if (ET_Group.find(group_id.substr(1)) != ET_Group.end() ) {
                        if (ET_SplitGroup) {
                            if (!ET_SplitType) {
                                ET_GroupStream.open((ET_OutputPrefix + group_id + "_extracted_data.fa").c_str());
                            } else {
                                if (ET_Spacer) {
                                    ET_SpacerStream.open((ET_OutputPrefix + group_id + "_spacers.fa").c_str());
                                }
                                if (ET_DirectRepeat) {
                                    ET_RepeatStream.open((ET_OutputPrefix + group_id + "_direct_repeats.fa").c_str());
                                }
                                if (ET_Flanker) {
                                    ET_FlankerStream.open((ET_OutputPrefix + group_id + "_flankers.fa").c_str());
                                }
                            }
                        }
                        
                        // matches to one of our wanted groups
                        extractDataFromGroup(xmlObj, currentElement);
                    } 
                    xr(&c_group_id);
                } else {
                    if (ET_SplitGroup) {
                        if (!ET_SplitType) {
                            ET_GroupStream.open((ET_OutputPrefix + group_id + "_extracted_data.fa").c_str());
                        } else {
                            if (ET_Spacer) {
                                ET_SpacerStream.open((ET_OutputPrefix + group_id + "_spacers.fa").c_str());
                            }
                            if (ET_DirectRepeat) {
                                ET_RepeatStream.open((ET_OutputPrefix + group_id + "_direct_repeats.fa").c_str());
                            }
                            if (ET_Flanker) {
                                ET_FlankerStream.open((ET_OutputPrefix + group_id + "_flankers.fa").c_str());
                            }
                        }
                    }

                    extractDataFromGroup(xmlObj, currentElement);
                    
                    if(ET_SplitGroup) {
                        // we are only spliting on type make three generic files
                        if (ET_SplitType) {
                            if (ET_Spacer) {
                                ET_SpacerStream.close();
                            }
                            if (ET_DirectRepeat) {
                                ET_RepeatStream.close();
                            }
                            if (ET_Flanker) {
                                ET_FlankerStream.close();
                            }
                        } else {
                            ET_GroupStream.close();
                        }

                    }
                }
            }
            
        }
        
    } catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        throw (crispr::xml_exception(__FILE__, __LINE__,__PRETTY_FUNCTION__,(errBuf.str()).c_str()));
        xercesc::XMLString::release( &message );
        
    } catch (crispr::xml_exception& xe) {
        std::cerr<< xe.what()<<std::endl;
        return;
    }
}

void ExtractTool::extractDataFromGroup(crispr::XML& xmlDoc, xercesc::DOMElement * currentGroup)
{
    // get the first child - the data element
    try {
        for (xercesc::DOMElement * currentElement = currentGroup->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
            char * c_gid = tc(currentGroup->getAttribute(xmlDoc.getGid()));
            if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc.getDrs())) {
                if (ET_DirectRepeat) {
                    // get direct repeats
                    if (ET_SplitType) {
                        processData(xmlDoc, currentElement, REPEAT, c_gid, ET_RepeatStream);
                    } else if (ET_SplitGroup) {
                        processData(xmlDoc, currentElement, REPEAT, c_gid, ET_GroupStream);   
                    } else {
                        processData(xmlDoc, currentElement, REPEAT, c_gid, ET_OneStream);   
                    }
                }
            } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc.getSpacers())) {
                if (ET_Spacer) {
                    // get spacers
                    if (ET_SplitType) {
                        processData(xmlDoc, currentElement, SPACER, c_gid, ET_SpacerStream);
                    } else if (ET_SplitGroup) {
                        processData(xmlDoc, currentElement, SPACER, c_gid, ET_GroupStream);   
                    } else {
                        processData(xmlDoc, currentElement, SPACER, c_gid, ET_OneStream);
                    }
                }
            } else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc.getFflankers())) {
                if (ET_Flanker) {
                    // get flankers
                    if (ET_SplitType) {
                        processData(xmlDoc, currentElement, FLANKER, c_gid, ET_FlankerStream);
                    } else if (ET_SplitGroup) {
                        processData(xmlDoc, currentElement, FLANKER, c_gid, ET_GroupStream);
                    } else {
                        processData(xmlDoc, currentElement, FLANKER, c_gid, ET_OneStream);

                    }
                }
            }
            xr(&c_gid);
            
        }
    } catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        xercesc::XMLString::release( &message );
    } catch (crispr::xml_exception& xe) {
        std::cerr<< xe.what()<<std::endl;
        return;
    }
}

void ExtractTool::processData(crispr::XML& xmlDoc, xercesc::DOMElement * currentType, ELEMENT_TYPE wantedType, std::string gid, std::ostream& outStream)
{
    try {
        for (xercesc::DOMElement * currentElement = currentType->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) {
            
            char * c_seq = tc(currentElement->getAttribute(xmlDoc.getSeq()));
            std::string id;
            
            switch (wantedType) {
                case REPEAT:
                {
                    char * c_id = tc(currentElement->getAttribute(xmlDoc.getDrid()));
                    id = c_id;
                    xr(&c_id);
                    break;
                }
                case SPACER:
                {
                    char * c_id = tc(currentElement->getAttribute(xmlDoc.getSpid()));
                    id = c_id;
                    if (currentElement->hasAttribute(xmlDoc.getCov())) {
                        char * c_cov = tc(currentElement->getAttribute(xmlDoc.getCov()));
                        id += "_Cov_"; 
                        id += c_cov;
                        xr(&c_cov);
                    }
                    xr(&c_id);
                    break;
                }
                case FLANKER:
                {
                    char * c_id = tc(currentElement->getAttribute(xmlDoc.getFlid()));
                    id = c_id;
                    xr(&c_id);
                    break;
                }
                case CONSENSUS:
                {
                    break;
                }
                default:
                {
                    throw (crispr::runtime_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Input element enum unknown"));
                    break;
                }
            }
            
            outStream<<'>'<<gid<<id<<std::endl<<c_seq<<std::endl;
            xr(&c_seq);
            
        }
        
    } catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        throw (crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, (errBuf.str()).c_str()));
        xercesc::XMLString::release( &message );
        
    } catch (crispr::xml_exception& xe) {
        std::cerr<< xe.what()<<std::endl;
        return;
        
    } catch (crispr::runtime_exception& re) {
        std::cerr<<re.what()<<std::endl;
        return;
        
    }
}

int extractMain (int argc, char ** argv)
{
    try {
		ExtractTool et;
		int opt_index = et.processOptions (argc, argv);
		if (opt_index >= argc) {
			throw crispr::input_exception("No input file provided" );
		
        } else {
			// get cracking and process that file
			return et.processInputFile(argv[opt_index]);
		}
	} catch(crispr::input_exception& re) {
        std::cerr<<re.what()<<std::endl;
        extractUsage();
        return 1;
    } catch(crispr::exception& ce ) {
		std::cerr<<ce.what()<<std::endl;
		return 1;
	}

}


void extractUsage (void)
{
	std::cout<<CRISPRTOOLS_PACKAGE_NAME<<" extract [-ghyxsdf] file.crispr"<<std::endl;
	std::cout<<"Options:"<<std::endl;
	std::cout<<"-h					print this handy help message"<<std::endl;
    std::cout<<"-o DIR              output file directory  [default: .]" <<std::endl; 
	std::cout<<"-g INT[,n]          a comma separated list of group IDs that you would like to extract data from."<<std::endl;
	std::cout<<"					Note that only the group number is needed, do not use prefixes like 'Group' or 'G', which"<<std::endl;
	std::cout<<"					are sometimes used in file names or in a .crispr file"<<std::endl;
	std::cout<<"-s					Extract the spacers of the listed group"<<std::endl;
	std::cout<<"-d					Extract the direct repeats of the listed group"<<std::endl;
	std::cout<<"-f					Extract the flanking sequences of the listed group"<<std::endl;
	std::cout<<"-x					Split the results into different files for each group.  If multiple types are set i.e. -sd"<<std::endl;
	std::cout<<"					then both the spacers and direct repeats from each group will be in the one file"<<std::endl;
	std::cout<<"-y					Split the results into different files for each type of sequence from all selected groups."<<std::endl;
	std::cout<<"					Only has an effect if multiple types are set."<<std::endl;

}
				
				
