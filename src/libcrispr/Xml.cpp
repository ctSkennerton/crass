// File: XML.cpp
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
#include <map>
#include <set>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

// local includes
#include "Xml.h"
#include "Exception.h"
#include "StringCheck.h"
#include "StlExt.h"

crispr::xml::base::base()
{
    try  
    {
        // Initialize Xerces infrastructure   
        this->init();
        // transcode all the member variables
        this->alloc();
    }
    catch( xercesc::XMLException& e )
    {
      char * message = xercesc::XMLString::transcode( e.getMessage() );
      std::cerr << "XML toolkit initialization error: " << message << std::endl;
      xercesc::XMLString::release( &message );
      // throw exception here to return ERROR_XERCES_INIT
    }
}


crispr::xml::base::~base(void)
{
    try 
    {
        // Free memory
        this->dealloc();
    }
    catch( ... )
    {
        std::cerr << "Unknown exception encountered in TagNamesdtor" << std::endl;
    }

    xercesc::XMLPlatformUtils::Terminate();  // Terminate after release of memory
}
void crispr::xml::base::init(void) {
    xercesc::XMLPlatformUtils::Initialize();
}
void crispr::xml::base::alloc(void) {
    // USE sed
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%TAG_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%ATTR_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    TAG_assembly = xercesc::XMLString::transcode("assembly");
    TAG_bf = xercesc::XMLString::transcode("bf");
    TAG_bflankers = xercesc::XMLString::transcode("bflankers");
    TAG_bs = xercesc::XMLString::transcode("bs");
    TAG_bspacers = xercesc::XMLString::transcode("bspacers");
    TAG_command = xercesc::XMLString::transcode("command");
    TAG_consensus = xercesc::XMLString::transcode("consensus");
    TAG_contig = xercesc::XMLString::transcode("contig");
    TAG_crispr = xercesc::XMLString::transcode("crispr");
    TAG_cspacer = xercesc::XMLString::transcode("cspacer");
    TAG_data = xercesc::XMLString::transcode("data");
    TAG_dr = xercesc::XMLString::transcode("dr");
    TAG_drs = xercesc::XMLString::transcode("drs");
    TAG_epos = xercesc::XMLString::transcode("epos");
    TAG_ff = xercesc::XMLString::transcode("ff");
    TAG_fflankers = xercesc::XMLString::transcode("fflankers");
    TAG_file = xercesc::XMLString::transcode("file");
    TAG_flanker = xercesc::XMLString::transcode("flanker");
    TAG_flankers = xercesc::XMLString::transcode("flankers");
    TAG_fs = xercesc::XMLString::transcode("fs");
    TAG_fspacers = xercesc::XMLString::transcode("fspacers");
    TAG_group = xercesc::XMLString::transcode("group");
    TAG_metadata = xercesc::XMLString::transcode("metadata");
    TAG_name = xercesc::XMLString::transcode("name");
    TAG_notes = xercesc::XMLString::transcode("notes");
    TAG_program = xercesc::XMLString::transcode("program");
    TAG_source = xercesc::XMLString::transcode("source");
    TAG_sources = xercesc::XMLString::transcode("sources");
    TAG_spacer = xercesc::XMLString::transcode("spacer");
    TAG_spacers = xercesc::XMLString::transcode("spacers");
    TAG_spos = xercesc::XMLString::transcode("spos");
    TAG_version = xercesc::XMLString::transcode("version");
    
    ATTR_accession = xercesc::XMLString::transcode("accession");
    ATTR_cid = xercesc::XMLString::transcode("cid");
    ATTR_confcnt = xercesc::XMLString::transcode("confcnt");
    ATTR_cov = xercesc::XMLString::transcode("cov");
    ATTR_directjoin = xercesc::XMLString::transcode("directjoin");
    ATTR_drconf = xercesc::XMLString::transcode("drconf");
    ATTR_drid = xercesc::XMLString::transcode("drid");
    ATTR_drseq = xercesc::XMLString::transcode("drseq");
    ATTR_flid = xercesc::XMLString::transcode("flid");
    ATTR_gid = xercesc::XMLString::transcode("gid");
    ATTR_seq = xercesc::XMLString::transcode("seq");
    ATTR_soid = xercesc::XMLString::transcode("soid");
    ATTR_spid = xercesc::XMLString::transcode("spid");
    ATTR_totcnt = xercesc::XMLString::transcode("totcnt");
    ATTR_type = xercesc::XMLString::transcode("type");
    ATTR_url = xercesc::XMLString::transcode("url");
    ATTR_version = xercesc::XMLString::transcode("version");
}

void crispr::xml::base::dealloc(void) {
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLString::release( \&TAG_\1 );%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLString::release( \&ATTR_\1 );%" | sort | uniq
    
    xercesc::XMLString::release( &TAG_assembly );
    xercesc::XMLString::release( &TAG_bf );
    xercesc::XMLString::release( &TAG_bflankers );
    xercesc::XMLString::release( &TAG_bs );
    xercesc::XMLString::release( &TAG_bspacers );
    xercesc::XMLString::release( &TAG_command );
    xercesc::XMLString::release( &TAG_consensus );
    xercesc::XMLString::release( &TAG_contig );
    xercesc::XMLString::release( &TAG_crispr );
    xercesc::XMLString::release( &TAG_cspacer );
    xercesc::XMLString::release( &TAG_data );
    xercesc::XMLString::release( &TAG_dr );
    xercesc::XMLString::release( &TAG_drs );
    xercesc::XMLString::release( &TAG_epos );
    xercesc::XMLString::release( &TAG_ff );
    xercesc::XMLString::release( &TAG_fflankers );
    xercesc::XMLString::release( &TAG_file );
    xercesc::XMLString::release( &TAG_flanker );
    xercesc::XMLString::release( &TAG_flankers );
    xercesc::XMLString::release( &TAG_fs );
    xercesc::XMLString::release( &TAG_fspacers );
    xercesc::XMLString::release( &TAG_group );
    xercesc::XMLString::release( &TAG_metadata );
    xercesc::XMLString::release( &TAG_name );
    xercesc::XMLString::release( &TAG_notes );
    xercesc::XMLString::release( &TAG_program );
    xercesc::XMLString::release( &TAG_source );
    xercesc::XMLString::release( &TAG_sources );
    xercesc::XMLString::release( &TAG_spacer );
    xercesc::XMLString::release( &TAG_spacers );
    xercesc::XMLString::release( &TAG_spos );
    xercesc::XMLString::release( &TAG_version );
    
    xercesc::XMLString::release( &ATTR_accession );
    xercesc::XMLString::release( &ATTR_cid );
    xercesc::XMLString::release( &ATTR_confcnt );
    xercesc::XMLString::release( &ATTR_cov );
    xercesc::XMLString::release( &ATTR_directjoin );
    xercesc::XMLString::release( &ATTR_drconf );
    xercesc::XMLString::release( &ATTR_drid );
    xercesc::XMLString::release( &ATTR_drseq );
    xercesc::XMLString::release( &ATTR_flid );
    xercesc::XMLString::release( &ATTR_gid );
    xercesc::XMLString::release( &ATTR_seq );
    xercesc::XMLString::release( &ATTR_soid );
    xercesc::XMLString::release( &ATTR_spid );
    xercesc::XMLString::release( &ATTR_totcnt );
    xercesc::XMLString::release( &ATTR_type );
    xercesc::XMLString::release( &ATTR_url );
    xercesc::XMLString::release( &ATTR_version );
}


crispr::xml::reader::reader()
{
    init();
    alloc();
    CX_FileParser = new xercesc::XercesDOMParser;
}

crispr::xml::reader::~reader()
{
    delete CX_FileParser;
    dealloc();
    
}
void crispr::xml::reader::parseXMLFile(std::string XMLFile, std::string& wantedGroup, std::string * directRepeat, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly)
{
    //-----
    // why not!
    //
    
    //logInfo("Parsing from " << XMLFile, 1);



    
    // Configure DOM parser.
    CX_FileParser->setValidationScheme( xercesc::XercesDOMParser::Val_Never );
    CX_FileParser->setDoNamespaces( false );
    CX_FileParser->setDoSchema( false );
    CX_FileParser->setLoadExternalDTD( false );
    
    try
    {
        CX_FileParser->parse( XMLFile.c_str() );
        
        // no need to free this pointer - owned by the parent parser object
        xercesc::DOMDocument * xmlDoc = CX_FileParser->getDocument();
        
        // Get the top-level element: 
        xercesc::DOMElement * elementRoot = xmlDoc->getDocumentElement();
        if( !elementRoot ) throw(std::runtime_error( "empty XML document" ));
        
        // find our wanted group
        xercesc::DOMElement * wanted_group_element = getWantedGroupFromRoot(elementRoot, wantedGroup, directRepeat);
        if (!wanted_group_element) 
        {
            throw (std::runtime_error("Could not find the input group."));
        }
        
        // get the assembly node
        xercesc::DOMElement * assembly_element = parseGroupForAssembly(wanted_group_element);
        if (assembly_element == NULL) {
            throw (std::runtime_error("no assembly tag for group."));
        }
        // get the contigs and spacers
        parseAssemblyForContigIds(assembly_element, wantedContigs, spacersForAssembly);
        
    }
    catch( xercesc::XMLException& e )
    {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        xercesc::XMLString::release( &message );
    }
}

           
xercesc::DOMElement * crispr::xml::reader::getWantedGroupFromRoot(xercesc::DOMElement * parentNode, std::string& wantedGroup, std::string * directRepeat)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
        if (xercesc::XMLString::equals(currentElement->getTagName(), tag_Group()))
        {
            // new group
            // test if it's one that we want
            //std::cout << "Group_" << XMLCH_2_STR(element->getAttribute(ATTR_gid)) << "_" << XMLCH_2_STR(element->getAttribute(ATTR_drseq)) << ".fa" << std::endl;
            char * c_group_name = tc(currentElement->getAttribute(attr_Gid()));
            std::string current_group_name = c_group_name;
            xr(&c_group_name);
            if (current_group_name == wantedGroup) 
            {
                // get the length of the direct repeat
                char * c_dr = tc(currentElement->getAttribute(attr_Drseq()));
                *directRepeat = c_dr;
                return currentElement;
            }
        }
        
    }
    
    // we should theoretically never get here but if the xml is bad then it might just happen
    // or if the user has put in a group that doesn't exist by mistake
    return NULL;
}

xercesc::DOMElement * crispr::xml::reader::parseGroupForAssembly(xercesc::DOMElement* parentNode)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
       if( xercesc::XMLString::equals(currentElement->getTagName(), tag_Assembly()))
       {
           // assembly section
           // the child nodes will be the contigs
           return currentElement;
       }       
   }
    // if there is no assembly for this group
    return NULL;
} 

void crispr::xml::reader::parseAssemblyForContigIds(xercesc::DOMElement* parentNode, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
       if( xercesc::XMLString::equals(currentElement->getTagName(), tag_Contig()))
       {
           // check to see if the current contig is one that we want
           char * c_current_contig = tc(currentElement->getAttribute(attr_Cid()));
           std::string current_contig = c_current_contig;
           std::set<std::string>::iterator contig_iter = wantedContigs.find(current_contig);
           if( contig_iter != wantedContigs.end())
           {
               
               // get the spacers from the assembly
               getSpacerIdForAssembly(currentElement, spacersForAssembly);
           }
           xr(&c_current_contig);
       }
       
   }
}

void crispr::xml::reader::getSpacerIdForAssembly(xercesc::DOMElement* parentNode, std::list<std::string>& spacersForAssembly)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
        if( xercesc::XMLString::equals(currentElement->getTagName(), tag_Cspacer()))
        {
            char * c_cspacer = tc(currentElement->getAttribute(attr_Spid()));
            std::string str = c_cspacer;
            spacersForAssembly.push_back(str);
            xr(&c_cspacer);
        }
    }
}

//DOMDocument * crispr::xml::base::setFileParser(const char * XMLFile)
//{
//    // Configure DOM parser.
//    xercesc::XercesDOMParser * d->setValidationScheme( XercesDOMParser::Val_Never );
//    xercesc::XercesDOMParser * d->setDoNamespaces( false );
//    xercesc::XercesDOMParser * d->setDoSchema( false );
//    xercesc::XercesDOMParser * d->setLoadExternalDTD( false );
//    
//    try
//    {
//        CX_FileParser->parse( XMLFile );
//        return CX_FileParser->getDocument();        
//    }
//    catch( xercesc::XMLException& e ) {
//        char* message = xercesc::XMLString::transcode( e.getMessage() );
//        std::stringstream errBuf;
//        errBuf << "Error parsing file: " << message << std::flush;
//        xercesc::XMLString::release( &message );
//        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
//    } catch (xercesc::DOMException& e) {
//        char* message = xercesc::XMLString::transcode( e.getMessage() );
//        std::stringstream errBuf;
//        errBuf << "Error parsing file: " << message << std::flush;
//        xercesc::XMLString::release( &message );
//        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
//    }
//}

crispr::xml::writer::writer() {
    XW_DocElem = NULL;
    init();
    alloc();
    
}

crispr::xml::writer::~writer() {
    dealloc();
}

xercesc::DOMElement * crispr::xml::writer::createDOMDocument(std::string rootElement, std::string versionNumber, int& errorNumber )   
{
    XMLCh * core = tc("Core");
    xercesc::DOMImplementation* impl =  xercesc::DOMImplementationRegistry::getDOMImplementation(core);
    xr(&core);
    if (impl != NULL)
    {
        try
        {
            XMLCh * x_root_elem = tc(rootElement.c_str());
            XW_DocElem = impl->createDocument( 0, x_root_elem, 0);  
            xr(&x_root_elem);

            if (XW_DocElem != NULL) 
            {
                xercesc::DOMElement* rootElem = XW_DocElem->getDocumentElement();
                XMLCh * x_version_num = tc(versionNumber.c_str());

                rootElem->setAttribute(attr_Version(), x_version_num);
                xr(&x_version_num);

                errorNumber = 0;
                return rootElem;
            }
        }
        catch (const xercesc::OutOfMemoryException&)
        {
            XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
            errorNumber =  5;
        }
        catch (const xercesc::DOMException& e)
        {
            XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
            errorNumber =  2;
        }
        catch (...)
        {
            XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
            errorNumber =  3;
        }
    }  // (inpl != NULL)
    else
    {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorNumber =  4;
    }
    return NULL;
}

xercesc::DOMElement * crispr::xml::writer::createDOMDocument(const char * rootElement, const char * versionNumber, int& errorNumber )   
{
    XMLCh * core = tc("Core");
    xercesc::DOMImplementation* impl =  xercesc::DOMImplementationRegistry::getDOMImplementation(core);
    xr(&core);
    if (impl != NULL)
    {
        try
        {
            XMLCh * x_root_elem = tc(rootElement);
            XW_DocElem = impl->createDocument(0, x_root_elem, 0);
            
            xr(&x_root_elem);
            
            if (XW_DocElem != NULL) 
            {
                xercesc::DOMElement* rootElem = XW_DocElem->getDocumentElement();
                XMLCh * x_version_num = tc(versionNumber);
                
                rootElem->setAttribute(attr_Version(), x_version_num );
                
                xr(&x_version_num);
                errorNumber = 0;
                return rootElem;
            }
        }
        catch (const xercesc::OutOfMemoryException&)
        {
            XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
            errorNumber =  5;
        }
        catch (const xercesc::DOMException& e)
        {
            XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
            errorNumber =  2;
        }
        catch (...)
        {
            XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
            errorNumber =  3;
        }
    }  // (inpl != NULL)
    else
    {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorNumber =  4;
    }
    return NULL;
}

xercesc::DOMElement * crispr::xml::writer::addMetaData(xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * meta_data_elem = XW_DocElem->createElement(tag_Metadata());
    parentNode->appendChild(meta_data_elem);    
    return meta_data_elem;
    
}

xercesc::DOMElement * crispr::xml::writer::addMetaData(std::string notes, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * meta_data_elem = XW_DocElem->createElement(tag_Metadata());
    xercesc::DOMElement * notes_elem = XW_DocElem->createElement(tag_Notes());
    XMLCh * x_notes = tc(notes.c_str());
    xercesc::DOMText * meta_data_notes = XW_DocElem->createTextNode(x_notes);
    xr(&x_notes);
    notes_elem->appendChild(meta_data_notes);
    meta_data_elem->appendChild(notes_elem);
    parentNode->appendChild(meta_data_elem);
    
    return meta_data_elem;
    
}

void crispr::xml::writer::addFileToMetadata(std::string type, std::string url, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * file = XW_DocElem->createElement(tag_File());
    XMLCh * x_type = tc(type.c_str());
    XMLCh * x_url = tc(url.c_str());
    file->setAttribute(attr_Type(), x_type);
    file->setAttribute(attr_Url(), x_url);
    xr(&x_type);
    xr(&x_url);
    parentNode->appendChild(file);
}

void crispr::xml::writer::addNotesToMetadata(std::string notes, xercesc::DOMElement *parentNode)
{
    xercesc::DOMElement * notes_elem = XW_DocElem->createElement(tag_Notes());
    XMLCh * x_notes = tc(notes.c_str());
    xercesc::DOMText * meta_data_notes = XW_DocElem->createTextNode(x_notes);
    xr(&x_notes);
    notes_elem->appendChild(meta_data_notes);
    parentNode->appendChild(notes_elem);
}

xercesc::DOMElement * crispr::xml::writer::addGroup(std::string& gID, std::string& drConsensus, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * group = XW_DocElem->createElement(tag_Group());
    
    // Set the attributes of the group
    XMLCh * x_gID = tc(gID.c_str());
    XMLCh * x_drConsensus = tc(drConsensus.c_str());
    
    group->setAttribute(attr_Gid(), x_gID);
    group->setAttribute(attr_Drseq(), x_drConsensus);
    
    xr(&x_gID);
    xr(&x_drConsensus);
    // add the group to the parent (root element)
    parentNode->appendChild(group);
    return group;
}
xercesc::DOMElement * crispr::xml::writer::addData(xercesc::DOMElement * parentNode)
{
    // create the data node with spacer and drs as child elements
    xercesc::DOMElement * sources = XW_DocElem->createElement(tag_Sources());
    xercesc::DOMElement * data = XW_DocElem->createElement(tag_Data());
    xercesc::DOMElement * drs = XW_DocElem->createElement(tag_Drs());
    xercesc::DOMElement * spacers = XW_DocElem->createElement(tag_Spacers());
    data->appendChild(sources);
    data->appendChild(drs);
    data->appendChild(spacers);
    parentNode->appendChild(data);
    return data;
}
xercesc::DOMElement * crispr::xml::writer::addAssembly(xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * assembly = XW_DocElem->createElement(tag_Assembly());
    parentNode->appendChild(assembly);
    return assembly;
}
void crispr::xml::writer::addDirectRepeat(std::string& drid, std::string& seq, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * dr = XW_DocElem->createElement(tag_Dr());
    
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_drid = tc(drid.c_str());
    
    dr->setAttribute(attr_Seq(), x_seq);
    dr->setAttribute(attr_Drid(), x_drid);
    
    xr(&x_seq);
    xr(&x_drid);
    
    parentNode->appendChild(dr);
}
xercesc::DOMElement * crispr::xml::writer::addSpacer(std::string& seq, std::string& spid, xercesc::DOMElement * parentNode, std::string cov)
{
    xercesc::DOMElement * sp = XW_DocElem->createElement(tag_Spacer());
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_cov = tc(cov.c_str());
    XMLCh * x_spid = tc(spid.c_str());
    
    sp->setAttribute(attr_Seq(), x_seq);
    sp->setAttribute(attr_Spid(), x_spid);
    sp->setAttribute(attr_Cov(), x_cov);
    
    xr(&x_seq);
    xr(&x_cov);
    xr(&x_spid);
    
    parentNode->appendChild(sp);
    return sp;
}
xercesc::DOMElement * crispr::xml::writer::createFlankers(xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * flankers = XW_DocElem->createElement(tag_Flankers());
    parentNode->appendChild(flankers);
    return flankers;
}
xercesc::DOMElement * crispr::xml::writer::addFlanker(std::string& seq, std::string& flid, xercesc::DOMElement * parentNode)
{
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_flid = tc(flid.c_str());

    xercesc::DOMElement * flanker = XW_DocElem->createElement(tag_Flanker());
    flanker->setAttribute(attr_Seq(), x_seq);
    flanker->setAttribute(attr_Flid(), x_flid);
    
    xr(&x_seq);
    xr(&x_flid);
    
    parentNode->appendChild(flanker);
    return flanker;
}
xercesc::DOMElement * crispr::xml::writer::addContig(std::string& cid, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * contig = XW_DocElem->createElement(tag_Contig());
    
    XMLCh * x_cid = tc(cid.c_str());
    contig->setAttribute(attr_Cid(), x_cid);
    xr(&x_cid);
    parentNode->appendChild(contig);
    return contig;
}
void crispr::xml::writer::createConsensus(std::string& concensus, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * concensus_elem = XW_DocElem->createElement(tag_Consensus());
    XMLCh * x_consensus = tc(concensus.c_str());
    xercesc::DOMText * concensus_text = XW_DocElem->createTextNode(x_consensus);
    
    xr(&x_consensus);
    
    concensus_elem->appendChild(concensus_text);
    parentNode->appendChild(concensus_elem);
}
xercesc::DOMElement * crispr::xml::writer::addSpacerToContig(std::string& spid, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * cspacer = XW_DocElem->createElement(tag_Cspacer());
    XMLCh * x_spid = tc(spid.c_str());
    cspacer->setAttribute(attr_Spid(), x_spid);
    xr(&x_spid);
    parentNode->appendChild(cspacer);
    return cspacer;
}
xercesc::DOMElement * crispr::xml::writer::createSpacers(std::string tag)
{
    XMLCh * x_tag = tc(tag.c_str());
    xercesc::DOMElement * spacers = XW_DocElem->createElement(x_tag);
    xr(&x_tag);
    return spacers;
}

xercesc::DOMElement * crispr::xml::writer::createFlankers(std::string tag)
{
    return createSpacers(tag);
}

void crispr::xml::writer::addSpacer(std::string tag, std::string& spid, std::string& drid, std::string& drconf, xercesc::DOMElement * parentNode)
{
    XMLCh * x_tag = tc(tag.c_str());
    XMLCh * x_drid = tc(drid.c_str());
    XMLCh * x_drconf = tc(drconf.c_str());
    XMLCh * x_spid = tc(spid.c_str());

    xercesc::DOMElement * fs = XW_DocElem->createElement(x_tag);
    fs->setAttribute(attr_Drid(), x_drid);
    fs->setAttribute(attr_Drconf(), x_drconf);
    fs->setAttribute(attr_Spid(), x_spid);
    
    xr(&x_tag);
    xr(&x_drid);
    xr(&x_drconf);
    xr(&x_spid);
    
    parentNode->appendChild(fs);
}
void crispr::xml::writer::addFlanker(std::string tag, std::string& flid, std::string& drconf, std::string& directjoin, xercesc::DOMElement * parentNode)
{
    XMLCh * x_tag = tc(tag.c_str());
    XMLCh * x_flid = tc(flid.c_str());
    XMLCh * x_drconf = tc(drconf.c_str());
    XMLCh * x_directjoin = tc(directjoin.c_str());
    
    xercesc::DOMElement * bf = XW_DocElem->createElement(x_tag);
    bf->setAttribute(attr_Flid(), x_flid);
    bf->setAttribute(attr_Drconf(), x_drconf);
    bf->setAttribute(attr_Directjoin(), x_directjoin);
    
    xr(&x_tag);
    xr(&x_flid);
    xr(&x_drconf);
    xr(&x_directjoin);
    
    parentNode->appendChild(bf);
}

//create the sources tag for a group <sources>
xercesc::DOMElement * crispr::xml::writer::addSources(xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * sources = XW_DocElem->createElement(tag_Sources());
    parentNode->appendChild(sources);
    return sources;
}

// create a source tag for either the sources in <group> 
xercesc::DOMElement * crispr::xml::writer::addSource(std::string accession, std::string soid, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * source = XW_DocElem->createElement(tag_Source());
    XMLCh * x_acc = tc(accession.c_str());
    source->setAttribute(attr_Accession(), x_acc);
    XMLCh * x_soid = tc(soid.c_str());
    source->setAttribute(attr_Soid(), x_soid);
    parentNode->appendChild(source);
    xr(&x_acc);
    xr(&x_soid);
    return source;
                        
}

xercesc::DOMElement * crispr::xml::writer::addSpacerSource(std::string soid, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * source = XW_DocElem->createElement(tag_Source());
    XMLCh * x_soid = tc(soid.c_str());
    source->setAttribute(attr_Soid(), x_soid);
    parentNode->appendChild(source);
    xr(&x_soid);
    return source;
}


// add start and end positions for <source> in <spacer>
void crispr::xml::writer::addStartAndEndPos(std::string start, std::string end, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * start_tag = XW_DocElem->createElement(tag_Spos());
    xercesc::DOMElement * end_tag = XW_DocElem->createElement(tag_Epos());
    XMLCh * x_start_site = tc(start.c_str());
    XMLCh * x_end_site = tc(end.c_str());
    xercesc::DOMText * start_text = XW_DocElem->createTextNode(x_start_site);
    xercesc::DOMText * end_text = XW_DocElem->createTextNode(x_end_site);
    start_tag->appendChild(start_text);
    end_tag->appendChild(end_text);
    parentNode->appendChild(start_tag);
    parentNode->appendChild(end_tag);
    xr(&x_start_site);
    xr(&x_end_site);
}

// add a <program> tag to <metadata>
xercesc::DOMElement * crispr::xml::writer::addProgram(xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * program = XW_DocElem->createElement(tag_Program());
    parentNode->appendChild(program);
    return program;
}

//add a <name> tag to <program>
void crispr::xml::writer::addProgName(std::string progName, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * name_tag = XW_DocElem->createElement(tag_Name());
    XMLCh * x_prog_name = tc(progName.c_str());
    xercesc::DOMText * name_text = XW_DocElem->createTextNode(x_prog_name);
    name_tag->appendChild(name_text);
    parentNode->appendChild(name_tag);
    xr(&x_prog_name);
}

// add a <version> tag to <program>
void crispr::xml::writer::addProgVersion(std::string progVersion, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * version_tag = XW_DocElem->createElement(tag_Version());
    XMLCh * x_prog_version = tc(progVersion.c_str());
    xercesc::DOMText * version_text = XW_DocElem->createTextNode(x_prog_version);
    version_tag->appendChild(version_text);
    parentNode->appendChild(version_tag);
    xr(&x_prog_version);
}

//add a <command> tag to <program>
void crispr::xml::writer::addProgCommand(std::string progCommand, xercesc::DOMElement * parentNode)
{
    xercesc::DOMElement * command_tag = XW_DocElem->createElement(tag_Command());
    XMLCh * x_prog_version = tc(progCommand.c_str());
    xercesc::DOMText * command_text = XW_DocElem->createTextNode(x_prog_version);
    command_tag->appendChild(command_text);
    parentNode->appendChild(command_tag);
    xr(&x_prog_version);
}

//
 // File IO / Printing
//

bool crispr::xml::writer::printDOMToFile(std::string outFileName )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {xercesc::chLatin_L, xercesc::chLatin_S, xercesc::chNull};
        xercesc::DOMImplementation *impl          = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
        xercesc::DOMLSSerializer   *theSerializer = ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();
        xercesc::DOMLSOutput       *theOutputDesc = ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);

        xercesc::DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        xercesc::XMLFormatTarget *myFormTarget;
        myFormTarget = new xercesc::LocalFileFormatTarget(outFileName.c_str());
        //myFormTarget=new StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);

        theSerializer->write(XW_DocElem, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const xercesc::OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (xercesc::XMLException& e)
    {
        char * c_exept = tc(e.getMessage());
        XERCES_STD_QUALIFIER cerr << "An error occurred during creation of output transcoder. Msg is:"
        << XERCES_STD_QUALIFIER endl
        << c_exept << XERCES_STD_QUALIFIER endl;
        retval = false;
        xr(&c_exept);
    }

return retval;

}

bool crispr::xml::writer::printDOMToScreen(void )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {xercesc::chLatin_L, xercesc::chLatin_S, xercesc::chNull};
        xercesc::DOMImplementation *impl          = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
        xercesc::DOMLSSerializer   *theSerializer = ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();
        xercesc::DOMLSOutput       *theOutputDesc = ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        xercesc::DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        xercesc::XMLFormatTarget *myFormTarget;
        myFormTarget=new xercesc::StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);
        
        theSerializer->write(XW_DocElem, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const xercesc::OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (xercesc::XMLException& e)
    {
        char * c_exept = tc(e.getMessage());
        XERCES_STD_QUALIFIER cerr << "An error occurred during creation of output transcoder. Msg is:"
        << XERCES_STD_QUALIFIER endl
        << c_exept << XERCES_STD_QUALIFIER endl;
        retval = false;
        xr(&c_exept);
    }
    
    return retval;
    
}

bool crispr::xml::writer::printDOMToFile(std::string outFileName, xercesc::DOMDocument * docDOM )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {xercesc::chLatin_L, xercesc::chLatin_S, xercesc::chNull};
        xercesc::DOMImplementation *impl          = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
        xercesc::DOMLSSerializer   *theSerializer = ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();
        xercesc::DOMLSOutput       *theOutputDesc = ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        
        xercesc::DOMConfiguration* serializerConfig = theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        xercesc::XMLFormatTarget *myFormTarget;
        myFormTarget = new xercesc::LocalFileFormatTarget(outFileName.c_str());
        //myFormTarget=new StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);
        
        theSerializer->write(docDOM, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const xercesc::OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (xercesc::XMLException& e)
    {
        char * c_msg = tc(e.getMessage());
        XERCES_STD_QUALIFIER cerr << "An error occurred during creation of output transcoder. Msg is:"
        << XERCES_STD_QUALIFIER endl
        << c_msg << XERCES_STD_QUALIFIER endl;
        retval = false;
        xr(&c_msg);
    }
    
    return retval;
    
}

bool crispr::xml::writer::printDOMToScreen(xercesc::DOMDocument * domDoc )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {xercesc::chLatin_L, xercesc::chLatin_S, xercesc::chNull};
        xercesc::DOMImplementation *impl          = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
        xercesc::DOMLSSerializer   *theSerializer = ((xercesc::DOMImplementationLS*)impl)->createLSSerializer();
        xercesc::DOMLSOutput       *theOutputDesc = ((xercesc::DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");

        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        
        xercesc::DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(xercesc::XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(xercesc::XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        xercesc::XMLFormatTarget *myFormTarget;
        myFormTarget=new xercesc::StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);
        
        theSerializer->write(domDoc, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const xercesc::OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (xercesc::XMLException& e)
    {
        char * c_msg = tc(e.getMessage());
        XERCES_STD_QUALIFIER cerr << "An error occurred during creation of output transcoder. Msg is:"
        << XERCES_STD_QUALIFIER endl
        << c_msg << XERCES_STD_QUALIFIER endl;
        retval = false;
        xr(&c_msg);
    }
    
    return retval;
    
}


