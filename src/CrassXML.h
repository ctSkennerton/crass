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


#include "StlExt.h"


// Error codes
enum {
    ERROR_ARGS = 1,
    ERROR_XERCES_INIT,
    ERROR_PARSE,
    ERROR_EMPTY_DOCUMENT
};

// ---------------------------------------------------------------------------
//  This is a simple class that lets us do easy (though not terribly efficient)
//  trancoding of char* data to XMLCh data.
// ---------------------------------------------------------------------------
class TranscodeStr
{
    public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    TranscodeStr(const char* const toTranscode)
    {
        // Call the private transcoding method
        TR_XmlString = xercesc::XMLString::transcode(toTranscode);
    }
    TranscodeStr(const std::string& toTranscode)
    {
        // Call the private transcoding method
        TR_XmlString = xercesc::XMLString::transcode(toTranscode.c_str());
    }
    TranscodeStr(const XMLCh * toTranscode)
    {
        // Call the private transcoding method
        TR_CString = xercesc::XMLString::transcode(toTranscode);
        TR_CUsed = true;
    }
    ~TranscodeStr()
    {
        if (TR_CUsed) {
            xercesc::XMLString::release(&TR_CString);
        } else {
            xercesc::XMLString::release(&TR_XmlString);
        }
    }
    
    
    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const char * cForm() const
    {
        if (TR_CUsed) {
            return TR_CString;
        } else {
            return NULL;
            
        }
    }
    const XMLCh * xForm() const 
    {
        if (!TR_CUsed) {
            return TR_XmlString;
        } else {
            return NULL;
        }
    }
    
    private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fUnicodeForm
    //      This is the Unicode XMLCh format of the string.
    // -----------------------------------------------------------------------
    XMLCh *   TR_XmlString;
    char *    TR_CString;
    bool      TR_CUsed;
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
    // grep ATTLIST crass.dtd | perl -ne 's/[^ ]* [^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * get$_(void) { return ATTR_$original; };\n";' | sort | uniq
    // grep ELEMENT crass.dtd | perl -ne 's/[^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * get$_(void) { return TAG_$original; };\n";' |sort | uniq
    
    inline XMLCh * getCid(void){return ATTR_cid;};
    inline XMLCh * getConfcnt(void){return ATTR_confcnt;};
    inline XMLCh * getCov(void){return ATTR_cov;};
    inline XMLCh * getDirectjoin(void){return ATTR_directjoin;};
    inline XMLCh * getDrconf(void){return ATTR_drconf;};
    inline XMLCh * getDrid(void){return ATTR_drid;};
    inline XMLCh * getDrseq(void){return ATTR_drseq;};
    inline XMLCh * getFlid(void){return ATTR_flid;};
    inline XMLCh * getGid(void){return ATTR_gid;};
    inline XMLCh * getSeq(void){return ATTR_seq;};
    inline XMLCh * getSpid(void){return ATTR_spid;};
    inline XMLCh * getTotcnt(void){return ATTR_totcnt;};
    inline XMLCh * getType(void){return ATTR_type;};
    inline XMLCh * getUrl(void){return ATTR_url;};
    inline XMLCh * getVersion(void){return ATTR_version;};
    
    inline XMLCh * getAssembly(void){return TAG_assembly;};
    inline XMLCh * getBf(void){return TAG_bf;};
    inline XMLCh * getBflankers(void){return TAG_bflankers;};
    inline XMLCh * getBs(void){return TAG_bs;};
    inline XMLCh * getBspacers(void){return TAG_bspacers;};
    inline XMLCh * getConsensus(void){return TAG_consensus;};
    inline XMLCh * getContig(void){return TAG_contig;};
    inline XMLCh * getCrispr(void){return TAG_crispr;};
    inline XMLCh * getCspacer(void){return TAG_cspacer;};
    inline XMLCh * getData(void){return TAG_data;};
    inline XMLCh * getDr(void){return TAG_dr;};
    inline XMLCh * getDrs(void){return TAG_drs;};
    inline XMLCh * getFf(void){return TAG_ff;};
    inline XMLCh * getFflankers(void){return TAG_fflankers;};
    inline XMLCh * getFile(void){return TAG_file;};
    inline XMLCh * getFlanker(void){return TAG_flanker;};
    inline XMLCh * getFlankers(void){return TAG_flankers;};
    inline XMLCh * getFs(void){return TAG_fs;};
    inline XMLCh * getFspacers(void){return TAG_fspacers;};
    inline XMLCh * getGroup(void){return TAG_group;};
    inline XMLCh * getMetadata(void){return TAG_metadata;};
    inline XMLCh * getNotes(void){return TAG_notes;};
    inline XMLCh * getSpacer(void){return TAG_spacer;};
    inline XMLCh * getSpacers(void){return TAG_spacers;};
    
    xercesc::DOMElement * getRootElement(void)
    {
        return CX_DocElem->getDocumentElement();
    }
    
    xercesc::DOMDocument * getDocumentObj(void)
    {
        return CX_DocElem;
    }
    
    //
    // Working functions
    //
    void parseCrassXMLFile(std::string XMLFile);
    void parseCrassXMLFile(std::string XMLFile, std::string& wantedGroup, std::string * directRepeat, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly);
    //std::string XMLCH_2_STR(const XMLCh* xmlch);
    
    char * XMLCH_2_STR( const XMLCh* toTranscode ) 
    {
        //-----
        // transcode on the fly, make sure to keep the values so we can release
        //
        char* ret_val = xercesc::XMLString::transcode(toTranscode); 
        CHAR_transcodes.push_back(ret_val);
        return ret_val;
    }
    
    XMLCh * STR_2_XMLCH( const std::string& toTranscode ) 
    {  
        //-----
        // transcode on the fly, make sure to keep the values so we can release
        //
        XMLCh * ret_val = xercesc::XMLString::transcode(toTranscode.c_str()); 
        XML_transcodes.push_back(ret_val);
        return ret_val;
    }
    xercesc::DOMElement * getWantedGroupFromRoot(xercesc::DOMElement * currentElement, std::string& wantedGroup, std::string * directRepeat);
    xercesc::DOMElement * parseGroupForAssembly(xercesc::DOMElement* currentElement);
    void parseAssemblyForContigIds(xercesc::DOMElement* currentElement, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly);
    void getSpacerIdForAssembly(xercesc::DOMElement* currentElement, std::list<std::string>& spacersForAssembly);
    
    
    // Parsing functions
    xercesc::DOMDocument * setFileParser(const char * XMLFile);
    
    
    
    //DOMDocument Creation returns root node
    xercesc::DOMElement * createDOMDocument(std::string& rootElement, std::string& versionNumber, int& errorNumber );
    xercesc::DOMElement * createDOMDocument(const char * rootElement, const char * versionNumber, int& errorNumber );
    
    // add a <metadata> tag to <group> with notes
    xercesc::DOMElement * addMetaData(std::string notes, xercesc::DOMElement * parentNode);
    
    // add in extra files to this group such as a reads file
    void addFileToMetadata(std::string type, std::string url, xercesc::DOMElement * parentNode);
    
    // add a <group> to <crass_assem>
    xercesc::DOMElement * addGroup(std::string& gID, std::string& drConsensus, xercesc::DOMElement * parentNode);
    
    // add the <data> to <group>
    xercesc::DOMElement * addData(xercesc::DOMElement * parentNode);
    
    // add the <assembly> to <group>
    xercesc::DOMElement * addAssembly(xercesc::DOMElement * parentNode);
    
    // add a direct repeat to the <data> (<dr>)
    void addDirectRepeat(std::string& drid, std::string& seq, xercesc::DOMElement * dataNode);
    
    // add a spcaer to the <data> (<spacer>)
    void addSpacer(std::string& seq, std::string& spid, xercesc::DOMElement * dataNode);
    
    // create a <flankers> tag in <data>
    xercesc::DOMElement * createFlankers(xercesc::DOMElement * parentNode);
    
    // add a <flanker> to the <data> section (<flanker>)
    void addFlanker(std::string& seq, std::string& flid, xercesc::DOMElement * dataNode);
    
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
    
    
    //
    // File IO / printing
    //
    bool printDOMToFile(std::string outFileName );
    bool printDOMToScreen(void);
    bool printDOMToFile(std::string outFileName, xercesc::DOMDocument * domDoc );
    bool printDOMToScreen( xercesc::DOMDocument * domDoc);
    
private:
    xercesc::XercesDOMParser * CX_FileParser;			// parsing object
    xercesc::DOMDocument * CX_DocElem;
    
    // we need to keep track of what it is we are transcoding so we can release at the end
    std::vector<XMLCh*> XML_transcodes;
    std::vector<char*> CHAR_transcodes;
    
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLCh\* ATTR_\1;%" | sort | uniq
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLCh\* TAG_\1;%" | sort | uniq
    XMLCh* ATTR_cid;
    XMLCh* ATTR_confcnt;
    XMLCh* ATTR_directjoin;
    XMLCh* ATTR_drconf;
    XMLCh* ATTR_drid;
    XMLCh* ATTR_drseq;
    XMLCh* ATTR_flid;
    XMLCh* ATTR_gid;
    XMLCh* ATTR_seq;
    XMLCh* ATTR_spid;
    XMLCh* ATTR_cov;
    XMLCh* ATTR_totcnt;
    XMLCh* ATTR_type;
    XMLCh* ATTR_url;
    XMLCh* ATTR_version;
    
    XMLCh* TAG_assembly;
    XMLCh* TAG_bf;
    XMLCh* TAG_bflankers;
    XMLCh* TAG_bs;
    XMLCh* TAG_bspacers;
    XMLCh* TAG_consensus;
    XMLCh* TAG_contig;
    XMLCh* TAG_crispr;
    XMLCh* TAG_cspacer;
    XMLCh* TAG_data;
    XMLCh* TAG_dr;
    XMLCh* TAG_drs;
    XMLCh* TAG_ff;
    XMLCh* TAG_fflankers;
    XMLCh* TAG_file;
    XMLCh* TAG_flanker;
    XMLCh* TAG_flankers;
    XMLCh* TAG_fs;
    XMLCh* TAG_fspacers;
    XMLCh* TAG_group;
    XMLCh* TAG_log;
    XMLCh* TAG_metadata;
    XMLCh* TAG_notes;
    XMLCh* TAG_spacer;
    XMLCh* TAG_spacers;
    
};


#endif //CrassXML_h

