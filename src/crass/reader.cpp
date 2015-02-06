/*
 *  crass.h is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011 Connor Skennerton & Michael Imelfort. All rights reserved. 
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

#include "reader.h"

crispr::xml::reader::reader()
{
    XR_FileParser = new xercesc::XercesDOMParser;
}

crispr::xml::reader::~reader()
{
    delete XR_FileParser;
}


xercesc::DOMDocument * crispr::xml::reader::setFileParser(const char * XMLFile)
{
    // Configure DOM parser.
    XR_FileParser->setValidationScheme( xercesc::XercesDOMParser::Val_Never );
    XR_FileParser->setDoNamespaces( false );
    XR_FileParser->setDoSchema( false );
    XR_FileParser->setLoadExternalDTD( false );
    
    try
    {
        XR_FileParser->parse( XMLFile );
        return XR_FileParser->getDocument();        
    }
    catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        xercesc::XMLString::release( &message );
        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
    } catch (xercesc::DOMException& e) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        xercesc::XMLString::release( &message );
        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
    }
}
