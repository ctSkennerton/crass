// File: XML.h
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

#ifndef XML_h
#define XML_h

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
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include <string>
#include <stdexcept>
#include <set>
#include <list>
#include <vector>

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif
#include "StringCheck.h"
#define tc(buf) xercesc::XMLString::transcode(buf)
#define xr(buf) xercesc::XMLString::release(buf)

// Error codes
enum {
    ERROR_ARGS = 1,
    ERROR_XERCES_INIT,
    ERROR_PARSE,
    ERROR_EMPTY_DOCUMENT
};
namespace crispr {
    namespace xml {
        class base {
        public:
            //constructor / destructor
            base(void);
            ~base(void);
            
            // funtions used to keep constructor and destructor clean
            // and so the inherated classes can use them
            void init(void);
            void alloc(void);
            void dealloc(void);
            
            //
            // Generic get
            //
            // grep ATTLIST crass.dtd | perl -ne 's/[^ ]* [^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * attr_$_(void) { return ATTR_$original; }\n";' | sort | uniq
            // grep ELEMENT crass.dtd | perl -ne 's/[^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * tag_$_(void) { return TAG_$original; }\n";' |sort | uniq
            
            inline XMLCh * attr_Accession(void) { return ATTR_accession; }
            inline XMLCh * attr_Cid(void) { return ATTR_cid; }
            inline XMLCh * attr_Confcnt(void) { return ATTR_confcnt; }
            inline XMLCh * attr_Cov(void) { return ATTR_cov; }
            inline XMLCh * attr_Directjoin(void) { return ATTR_directjoin; }
            inline XMLCh * attr_Drconf(void) { return ATTR_drconf; }
            inline XMLCh * attr_Drid(void) { return ATTR_drid; }
            inline XMLCh * attr_Drseq(void) { return ATTR_drseq; }
            inline XMLCh * attr_Flid(void) { return ATTR_flid; }
            inline XMLCh * attr_Gid(void) { return ATTR_gid; }
            inline XMLCh * attr_Seq(void) { return ATTR_seq; }
            inline XMLCh * attr_Soid(void) { return ATTR_soid; }
            inline XMLCh * attr_Spid(void) { return ATTR_spid; }
            inline XMLCh * attr_Totcnt(void) { return ATTR_totcnt; }
            inline XMLCh * attr_Type(void) { return ATTR_type; }
            inline XMLCh * attr_Url(void) { return ATTR_url; }
            inline XMLCh * attr_Version(void) { return ATTR_version; }
            
            inline XMLCh * tag_Assembly(void) { return TAG_assembly; }
            inline XMLCh * tag_Bf(void) { return TAG_bf; }
            inline XMLCh * tag_Bflankers(void) { return TAG_bflankers; }
            inline XMLCh * tag_Bs(void) { return TAG_bs; }
            inline XMLCh * tag_Bspacers(void) { return TAG_bspacers; }
            inline XMLCh * tag_Command(void) { return TAG_command; }
            inline XMLCh * tag_Consensus(void) { return TAG_consensus; }
            inline XMLCh * tag_Contig(void) { return TAG_contig; }
            inline XMLCh * tag_Crispr(void) { return TAG_crispr; }
            inline XMLCh * tag_Cspacer(void) { return TAG_cspacer; }
            inline XMLCh * tag_Data(void) { return TAG_data; }
            inline XMLCh * tag_Dr(void) { return TAG_dr; }
            inline XMLCh * tag_Drs(void) { return TAG_drs; }
            inline XMLCh * tag_Epos(void) { return TAG_epos; }
            inline XMLCh * tag_Ff(void) { return TAG_ff; }
            inline XMLCh * tag_Fflankers(void) { return TAG_fflankers; }
            inline XMLCh * tag_File(void) { return TAG_file; }
            inline XMLCh * tag_Flanker(void) { return TAG_flanker; }
            inline XMLCh * tag_Flankers(void) { return TAG_flankers; }
            inline XMLCh * tag_Fs(void) { return TAG_fs; }
            inline XMLCh * tag_Fspacers(void) { return TAG_fspacers; }
            inline XMLCh * tag_Group(void) { return TAG_group; }
            inline XMLCh * tag_Metadata(void) { return TAG_metadata; }
            inline XMLCh * tag_Name(void) { return TAG_name; }
            inline XMLCh * tag_Notes(void) { return TAG_notes; }
            inline XMLCh * tag_Program(void) { return TAG_program; }
            inline XMLCh * tag_Source(void) { return TAG_source; }
            inline XMLCh * tag_Sources(void) { return TAG_sources; }
            inline XMLCh * tag_Spos(void) { return TAG_spos; }
            inline XMLCh * tag_Spacer(void) { return TAG_spacer; }
            inline XMLCh * tag_Spacers(void) { return TAG_spacers; }
            inline XMLCh * tag_Version(void) { return TAG_version; }
        
            
                        // Parsing functions
            xercesc::DOMDocument * setFileParser(const char * xmlFile);

        private:            
            // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLCh\* ATTR_\1;%" | sort | uniq
            // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLCh\* TAG_\1;%" | sort | uniq
            XMLCh * ATTR_accession;
            XMLCh * ATTR_cid;
            XMLCh * ATTR_confcnt;
            XMLCh * ATTR_cov;
            XMLCh * ATTR_directjoin;
            XMLCh * ATTR_drconf;
            XMLCh * ATTR_drid;
            XMLCh * ATTR_drseq;
            XMLCh * ATTR_flid;
            XMLCh * ATTR_gid;
            XMLCh * ATTR_seq;
            XMLCh * ATTR_soid;
            XMLCh * ATTR_spid;
            XMLCh * ATTR_totcnt;
            XMLCh * ATTR_type;
            XMLCh * ATTR_url;
            XMLCh * ATTR_version;
            
            XMLCh * TAG_assembly;
            XMLCh * TAG_bf;
            XMLCh * TAG_bflankers;
            XMLCh * TAG_bs;
            XMLCh * TAG_bspacers;
            XMLCh * TAG_command;
            XMLCh * TAG_consensus;
            XMLCh * TAG_contig;
            XMLCh * TAG_crispr;
            XMLCh * TAG_cspacer;
            XMLCh * TAG_data;
            XMLCh * TAG_dr;
            XMLCh * TAG_drs;
            XMLCh * TAG_epos;
            XMLCh * TAG_ff;
            XMLCh * TAG_fflankers;
            XMLCh * TAG_file;
            XMLCh * TAG_flanker;
            XMLCh * TAG_flankers;
            XMLCh * TAG_fs;
            XMLCh * TAG_fspacers;
            XMLCh * TAG_group;
            XMLCh * TAG_metadata;
            XMLCh * TAG_name;
            XMLCh * TAG_notes;
            XMLCh * TAG_program;
            XMLCh * TAG_source;
            XMLCh * TAG_sources;
            XMLCh * TAG_spacer;
            XMLCh * TAG_spacers;
            XMLCh * TAG_spos;
            XMLCh * TAG_version;
            
        };

        class reader : public base {
        public:
            reader();
            ~reader();
            
            //DOMDocument Creation returns root node
            void parseXMLFile(std::string XMLFile);
            //
            // Working functions
            //
            void parseXMLFile(std::string XMLFile, std::string& wantedGroup, std::string * directRepeat, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly);
            
            xercesc::DOMElement * getWantedGroupFromRoot(xercesc::DOMElement * currentElement, std::string& wantedGroup, std::string * directRepeat);
            xercesc::DOMElement * parseGroupForAssembly(xercesc::DOMElement* currentElement);
            void parseAssemblyForContigIds(xercesc::DOMElement* currentElement, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly);
            void getSpacerIdForAssembly(xercesc::DOMElement* currentElement, std::list<std::string>& spacersForAssembly);
            
            

            
        private:
            xercesc::XercesDOMParser * CX_FileParser;			// parsing object

        };
        
        class writer : public base {
            
            //members
            xercesc::DOMDocument * XW_DocElem;
            StringCheck XW_SourcesMap;
            int XW_CurrentSourceId;
            
        public:
            
            //constructor/destructor
            writer();
            ~writer();
            
            //
            // File IO / printing
            //
            xercesc::DOMElement * createDOMDocument(std::string rootElement, std::string versionNumber, int& errorNumber );
            
            xercesc::DOMElement * createDOMDocument(const char * rootElement, const char * versionNumber, int& errorNumber );
            
            // add a <metadata> tag to <group>
            xercesc::DOMElement * addMetaData(xercesc::DOMElement * parentNode);
            
            // add a <metadata> tag to <group> with notes
            xercesc::DOMElement * addMetaData(std::string notes, xercesc::DOMElement * parentNode);
            
            // add in extra files to this group such as a reads file
            void addFileToMetadata(std::string type, std::string url, xercesc::DOMElement * parentNode);
            
            void addNotesToMetadata(std::string notes, xercesc::DOMElement * parentNode);
            
            // add a <group> to <crass_assem>
            xercesc::DOMElement * addGroup(std::string& gID, std::string& drConsensus, xercesc::DOMElement * parentNode);
            
            // add the <data> to <group>
            xercesc::DOMElement * addData(xercesc::DOMElement * parentNode);
            
            // add the <assembly> to <group>
            xercesc::DOMElement * addAssembly(xercesc::DOMElement * parentNode);
            
            // add a direct repeat to the <data> (<dr>)
            void addDirectRepeat(std::string& drid, std::string& seq, xercesc::DOMElement * dataNode);
            
            // add a spcaer to the <data> (<spacer>)
            xercesc::DOMElement * addSpacer(std::string& seq, std::string& spid, xercesc::DOMElement * dataNode, std::string cov = "0" );
            
            // create a <flankers> tag in <data>
            xercesc::DOMElement * createFlankers(xercesc::DOMElement * parentNode);
            
            // add a <flanker> to the <data> section (<flanker>)
            xercesc::DOMElement * addFlanker(std::string& seq, std::string& flid, xercesc::DOMElement * dataNode);
            
            // add a <contig> to an <assembly>
            xercesc::DOMElement * addContig(std::string& cid, xercesc::DOMElement * parentNode);
            
            // add the concensus sequence to the contig (creates <concensus>)
            void createConsensus(std::string& concensus, xercesc::DOMElement * parentNode);
            
            // add a spacer to a contig (<cspacer> tag)
            xercesc::DOMElement * addSpacerToContig(std::string& spid, xercesc::DOMElement * parentNode);
            
            // use to create a backward spacers (<bspacers>) or forward spacers (<fspacers>) element by changing the value of tag
            xercesc::DOMElement * createSpacers(std::string tag);
            
            // use to create a backward flankers (<bflankers>) or forward flankers (<fflankers>) element by changing the value of tag
            xercesc::DOMElement * createFlankers(std::string tag);
            
            // use to add either a backward spacer (<bs>) or forward spacer (<fs>) by changing the value of tag
            void addSpacer(std::string tag, std::string& spid, std::string& drid, std::string& drconf, xercesc::DOMElement * parentNode);
            
            // use to add either a backward flanker (<bf>) or forward flanker (<ff>) by changing the value of tag
            void addFlanker(std::string tag, std::string& flid, std::string& drconf, std::string& directjoin, xercesc::DOMElement * parentNode);
            
            //create the sources tag for a group <sources>
            xercesc::DOMElement * addSources(xercesc::DOMElement * parentNode);
            
            // create a source tag for either the sources in <group> 
            xercesc::DOMElement * addSource(std::string accession, std::string soid, xercesc::DOMElement * parentNode);
            
            // create a source tag for the <spacer>
            xercesc::DOMElement * addSpacerSource(std::string soid, xercesc::DOMElement * parentNode);
            
            // add start and end positions for <source> in <spacer>
            void addStartAndEndPos(std::string start, std::string end, xercesc::DOMElement * parentNode);
            
            // add a <program> tag to <metadata>
            xercesc::DOMElement * addProgram(xercesc::DOMElement * parentNode);
            
            //add a <name> tag to <program>
            void addProgName(std::string progName, xercesc::DOMElement * parentNode);
            
            // add a <version> tag to <program>
            void addProgVersion(std::string progVersion, xercesc::DOMElement * parentNode);
            
            //add a <command> tag to <program>
            void addProgCommand(std::string progCommand, xercesc::DOMElement * parentNode);
            
            bool printDOMToFile(std::string outFileName );
            
            bool printDOMToScreen(void);
            
            bool printDOMToFile(std::string outFileName, xercesc::DOMDocument * domDoc );
            
            bool printDOMToScreen( xercesc::DOMDocument * domDoc);
            
            inline xercesc::DOMElement * getRootElement(void)
            {
                return XW_DocElem->getDocumentElement();
            }
            
            inline xercesc::DOMDocument * getDocumentObj(void)
            {
                return XW_DocElem;
            }
            
            
        };
    }
}


#endif //XML_h
