// File: CrassXML.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// Header file for the crass XML reader writer
//
// Many thanks to http://www.yolinux.com/TUTORIALS/XML-Xerces-C.html
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

#ifndef CrassXML_h
    #define CrassXML_h

// system includes
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/TransService.hpp>

#include <string>
#include <stdexcept>

// local includes
#include "crassDefines.h"

// Error codes
enum {
   ERROR_ARGS = 1,
   ERROR_XERCES_INIT,
   ERROR_PARSE,
   ERROR_EMPTY_DOCUMENT
};

class CrassXML 
{
    public:
        //constructor / destructor
        CrassXML(void);
        ~CrassXML(void);
        
        //
        // Generic get and set
        //
        
        //
        // Working functions
        //
        void parseCrassXMLFile(std::string XMLFile);
        std::string XMLCH_2_STR(const XMLCh* xmlch);

        //
        // File IO / printing
        //

    
    private:
        xercesc::XercesDOMParser * mConfigFileParser;			// parsing object
        // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLCh\* ATTR_\1;%" | sort | uniq
        // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLCh\* TAG_\1;%" | sort | uniq
        XMLCh* ATTR_cid;
        XMLCh* ATTR_confcnt;
        XMLCh* ATTR_drconf;
        XMLCh* ATTR_drid;
        XMLCh* ATTR_drseq;
        XMLCh* ATTR_flid;
        XMLCh* ATTR_gid;
        XMLCh* ATTR_seq;
        XMLCh* ATTR_spid;
        XMLCh* ATTR_totcnt;
        XMLCh* ATTR_version;
        XMLCh* TAG_assembly;
        XMLCh* TAG_bs;
        XMLCh* TAG_bspacers;
        XMLCh* TAG_contig;
        XMLCh* TAG_crass_assem;
        XMLCh* TAG_cspacer;
        XMLCh* TAG_data;
        XMLCh* TAG_dr;
        XMLCh* TAG_drs;
        XMLCh* TAG_flanker;
        XMLCh* TAG_flankers;
        XMLCh* TAG_fs;
        XMLCh* TAG_fspacers;
        XMLCh* TAG_group;
        XMLCh* TAG_log;
        XMLCh* TAG_spacer;
        XMLCh* TAG_spacers;
};

#endif //CrassXML_h
