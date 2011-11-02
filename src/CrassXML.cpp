// File: CrassXML.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
//  Crass XML reader writer. Write and read assemblies to and from file
//
// --------------------------------------------------------------------
//  Copyright  2011 Michael Imelfort and Connor Skennerton
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//
// system includes
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <list>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

// local includes
#include "CrassXML.h"
#include "LoggerSimp.h"

using namespace xercesc;

CrassXML::CrassXML(void)
{
    try  // Initialize Xerces infrastructure
    {
      XMLPlatformUtils::Initialize();
    }
    catch( XMLException& e )
    {
      char* message = XMLString::transcode( e.getMessage() );
      std::cerr << "XML toolkit initialization error: " << message << std::endl;
      XMLString::release( &message );
      // throw exception here to return ERROR_XERCES_INIT
    }
    mConfigFileParser = new XercesDOMParser;
    
    
    // gotta hard code all these motherfuckers
    // USE sed
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%TAG_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLString::release( \&TAG_\1 );%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%ATTR_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLString::release( \&ATTR_\1 );%" | sort | uniq
    TAG_assembly = XMLString::transcode("assembly");
    TAG_bspacers = XMLString::transcode("bspacers");
    TAG_bs = XMLString::transcode("bs");
    TAG_contig = XMLString::transcode("contig");
    TAG_crass_assem = XMLString::transcode("crass_assem");
    TAG_cspacer = XMLString::transcode("cspacer");
    TAG_data = XMLString::transcode("data");
    TAG_drs = XMLString::transcode("drs");
    TAG_dr = XMLString::transcode("dr");
    TAG_flankers = XMLString::transcode("flankers");
    TAG_flanker = XMLString::transcode("flanker");
    TAG_fspacers = XMLString::transcode("fspacers");
    TAG_fs = XMLString::transcode("fs");
    TAG_group = XMLString::transcode("group");
    TAG_log = XMLString::transcode("log");
    TAG_spacers = XMLString::transcode("spacers");
    TAG_spacer = XMLString::transcode("spacer");
    
    ATTR_cid = XMLString::transcode("cid");
    ATTR_confcnt = XMLString::transcode("confcnt");
    ATTR_drconf = XMLString::transcode("drconf");
    ATTR_drid = XMLString::transcode("drid");
    ATTR_drseq = XMLString::transcode("drseq");
    ATTR_flid = XMLString::transcode("flid");
    ATTR_gid = XMLString::transcode("gid");
    ATTR_seq = XMLString::transcode("seq");
    ATTR_spid = XMLString::transcode("spid");
    ATTR_totcnt = XMLString::transcode("totcnt");
    ATTR_version = XMLString::transcode("version");
}


CrassXML::~CrassXML(void)
{
    // Free memory
    delete mConfigFileParser;
    
    try // Free memory
    {
       XMLString::release( &TAG_assembly );
       XMLString::release( &TAG_bs );
       XMLString::release( &TAG_bspacers );
       XMLString::release( &TAG_contig );
       XMLString::release( &TAG_crass_assem );
       XMLString::release( &TAG_cspacer );
       XMLString::release( &TAG_data );
       XMLString::release( &TAG_dr );
       XMLString::release( &TAG_drs );
       XMLString::release( &TAG_flanker );
       XMLString::release( &TAG_flankers );
       XMLString::release( &TAG_fs );
       XMLString::release( &TAG_fspacers );
       XMLString::release( &TAG_group );
       XMLString::release( &TAG_log );
       XMLString::release( &TAG_spacer );
       XMLString::release( &TAG_spacers );
       
       XMLString::release( &ATTR_cid );
       XMLString::release( &ATTR_confcnt );
       XMLString::release( &ATTR_drconf );
       XMLString::release( &ATTR_drid );
       XMLString::release( &ATTR_drseq );
       XMLString::release( &ATTR_flid );
       XMLString::release( &ATTR_gid );
       XMLString::release( &ATTR_seq );
       XMLString::release( &ATTR_spid );
       XMLString::release( &ATTR_totcnt );
       XMLString::release( &ATTR_version );
    }
    catch( ... )
    {
        std::cerr << "Unknown exception encountered in TagNamesdtor" << std::endl;
    }
    
    try // Terminate Xerces
    {
      XMLPlatformUtils::Terminate();  // Terminate after release of memory
    }
    catch( xercesc::XMLException& e )
    {
      char* message = xercesc::XMLString::transcode( e.getMessage() );

      std::cerr << "XML ttolkit teardown error: " << message << std::endl;
      XMLString::release( &message );
    }
}

void CrassXML::parseCrassXMLFile(std::string XMLFile)
{
    //-----
    // why not!
    //

    logInfo("Parsing from " << XMLFile, 1);
    // Test to see if the file is ok.
    struct stat fileStatus;

    int iretStat = stat(XMLFile.c_str(), &fileStatus);
    if( iretStat == ENOENT )
       throw ( std::runtime_error("Path file_name does not exist, or path is an empty string.") );
    else if( iretStat == ENOTDIR )
       throw ( std::runtime_error("A component of the path is not a directory."));
    else if( iretStat == ELOOP )
       throw ( std::runtime_error("Too many symbolic links encountered while traversing the path."));
    else if( iretStat == EACCES )
       throw ( std::runtime_error("Permission denied."));
    else if( iretStat == ENAMETOOLONG )
       throw ( std::runtime_error("File can not be read\n"));

    // Configure DOM parser.
    mConfigFileParser->setValidationScheme( XercesDOMParser::Val_Never );
    mConfigFileParser->setDoNamespaces( false );
    mConfigFileParser->setDoSchema( false );
    mConfigFileParser->setLoadExternalDTD( false );

    try
       {
          mConfigFileParser->parse( XMLFile.c_str() );

          // no need to free this pointer - owned by the parent parser object
          DOMDocument* xmlDoc = mConfigFileParser->getDocument();

          // Get the top-level element: 
          DOMElement* elementRoot = xmlDoc->getDocumentElement();
          if( !elementRoot ) throw(std::runtime_error( "empty XML document" ));

          // get the children
          DOMNodeList*      children = elementRoot->getChildNodes();
          const  XMLSize_t nodeCount = children->getLength();

          // For all nodes, children of "root" in the XML tree.
          for( XMLSize_t xx = 0; xx < nodeCount; ++xx )
          {
             DOMNode* currentNode = children->item(xx);
             if( currentNode->getNodeType() &&  // true is not NULL
                 currentNode->getNodeType() == DOMNode::ELEMENT_NODE ) // is element 
             {
                // Found node which is an Element. Re-cast node as element
                DOMElement* currentElement
                            = dynamic_cast< xercesc::DOMElement* >( currentNode );
                if( XMLString::equals(currentElement->getTagName(), TAG_log))
                {
                    // log section
                }
                else if (XMLString::equals(currentElement->getTagName(), TAG_group))
                {
                    // new group
                    std::cout << XMLCH_2_STR(currentElement->getTagName()) << std::endl;
                    std::cout << "Group_" << XMLCH_2_STR(currentElement->getAttribute(ATTR_gid)) << "_" << XMLCH_2_STR(currentElement->getAttribute(ATTR_drseq)) << ".fa" << std::endl;
                    
                }
             }
          }
       }
       catch( xercesc::XMLException& e )
       {
          char* message = xercesc::XMLString::transcode( e.getMessage() );
          std::ostringstream errBuf;
          errBuf << "Error parsing file: " << message << std::flush;
          XMLString::release( &message );
       }
}

std::string CrassXML::XMLCH_2_STR(const XMLCh* xmlch)
{
    //-----
    // convert a XMCH* to a std::string in utf8
    //
    XMLTranscoder* utf8Transcoder;
    XMLTransService::Codes failReason;
    utf8Transcoder = XMLPlatformUtils::fgTransService->makeNewTranscoderFor("UTF-8", failReason, 16*1024);

    XMLSize_t len = XMLString::stringLen(xmlch);
    XMLByte* utf8 = new XMLByte(); // ?
    XMLSize_t eaten;
    XMLSize_t utf8Len = utf8Transcoder->transcodeTo(xmlch, len, utf8, len, eaten, XMLTranscoder::UnRep_Throw);

    utf8[utf8Len] = '\0';
    std::string str = (char*)utf8;
    
    delete utf8;
    delete utf8Transcoder;
    return str;
}
