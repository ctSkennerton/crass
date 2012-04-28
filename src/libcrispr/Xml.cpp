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

using namespace xercesc;

crispr::XML::XML(void)
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
    CX_FileParser = new XercesDOMParser;
    CX_DocElem = NULL;
    
    // USE sed
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%TAG_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLString::release( \&TAG_\1 );%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%ATTR_\1 = XMLString::transcode(\"\1\");%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLString::release( \&ATTR_\1 );%" | sort | uniq
    TAG_assembly = XMLString::transcode("assembly");
    TAG_bf = XMLString::transcode("bf");
    TAG_bflankers = XMLString::transcode("bflankers");
    TAG_bs = XMLString::transcode("bs");
    TAG_bspacers = XMLString::transcode("bspacers");
    TAG_consensus = XMLString::transcode("consensus");
    TAG_contig = XMLString::transcode("contig");
    TAG_crispr = XMLString::transcode("crispr");
    TAG_cspacer = XMLString::transcode("cspacer");
    TAG_data = XMLString::transcode("data");
    TAG_dr = XMLString::transcode("dr");
    TAG_drs = XMLString::transcode("drs");
    TAG_ff = XMLString::transcode("ff");
    TAG_fflankers = XMLString::transcode("fflankers");
    TAG_file = XMLString::transcode("file");
    TAG_flanker = XMLString::transcode("flanker");
    TAG_flankers = XMLString::transcode("flankers");
    TAG_fs = XMLString::transcode("fs");
    TAG_fspacers = XMLString::transcode("fspacers");
    TAG_group = XMLString::transcode("group");
    TAG_metadata = XMLString::transcode("metadata");
    TAG_notes = XMLString::transcode("notes");
    TAG_spacer = XMLString::transcode("spacer");
    TAG_spacers = XMLString::transcode("spacers");
    
    ATTR_cid = XMLString::transcode("cid");
    ATTR_confcnt = XMLString::transcode("confcnt");
    ATTR_directjoin = XMLString::transcode("directjoin");
    ATTR_drconf = XMLString::transcode("drconf");
    ATTR_drid = XMLString::transcode("drid");
    ATTR_drseq = XMLString::transcode("drseq");
    ATTR_flid = XMLString::transcode("flid");
    ATTR_gid = XMLString::transcode("gid");
    ATTR_seq = XMLString::transcode("seq");
    ATTR_spid = XMLString::transcode("spid");
    ATTR_cov = XMLString::transcode("cov");

    ATTR_totcnt = XMLString::transcode("totcnt");
    ATTR_type = XMLString::transcode("type");
    ATTR_url = XMLString::transcode("url");
    ATTR_version = XMLString::transcode("version");

}


crispr::XML::~XML(void)
{
    // Free memory
    delete CX_FileParser;
    if (CX_DocElem != NULL) 
    {
        CX_DocElem->release();
    }
    try // Free memory
    {
        
        XMLString::release( &TAG_assembly );
        XMLString::release( &TAG_bf );
        XMLString::release( &TAG_bflankers );
        XMLString::release( &TAG_bs );
        XMLString::release( &TAG_bspacers );
        XMLString::release( &TAG_consensus );
        XMLString::release( &TAG_contig );
        XMLString::release( &TAG_crispr );
        XMLString::release( &TAG_cspacer );
        XMLString::release( &TAG_data );
        XMLString::release( &TAG_dr );
        XMLString::release( &TAG_drs );
        XMLString::release( &TAG_ff );
        XMLString::release( &TAG_fflankers );
        XMLString::release( &TAG_file );
        XMLString::release( &TAG_flanker );
        XMLString::release( &TAG_flankers );
        XMLString::release( &TAG_fs );
        XMLString::release( &TAG_fspacers );
        XMLString::release( &TAG_group );
        XMLString::release( &TAG_metadata );
        XMLString::release( &TAG_notes );
        XMLString::release( &TAG_spacer );
        XMLString::release( &TAG_spacers );
        
        XMLString::release( &ATTR_cid );
        XMLString::release( &ATTR_confcnt );
        XMLString::release( &ATTR_directjoin );
        XMLString::release( &ATTR_drconf );
        XMLString::release( &ATTR_drid );
        XMLString::release( &ATTR_drseq );
        XMLString::release( &ATTR_flid );
        XMLString::release( &ATTR_gid );
        XMLString::release( &ATTR_seq );
        XMLString::release( &ATTR_spid );
        XMLString::release( &ATTR_cov );

        XMLString::release( &ATTR_totcnt );
        XMLString::release( &ATTR_type );
        XMLString::release( &ATTR_url );
        XMLString::release( &ATTR_version );
        

    }
    catch( ... )
    {
        std::cerr << "Unknown exception encountered in TagNamesdtor" << std::endl;
    }

    XMLPlatformUtils::Terminate();  // Terminate after release of memory
}


void crispr::XML::parseXMLFile(std::string XMLFile, std::string& wantedGroup, std::string * directRepeat, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly)
{
    //-----
    // why not!
    //
    
    //logInfo("Parsing from " << XMLFile, 1);



    
    // Configure DOM parser.
    CX_FileParser->setValidationScheme( XercesDOMParser::Val_Never );
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
        XMLString::release( &message );
    }
}

           
xercesc::DOMElement * crispr::XML::getWantedGroupFromRoot(xercesc::DOMElement * parentNode, std::string& wantedGroup, std::string * directRepeat)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
        if (XMLString::equals(currentElement->getTagName(), TAG_group))
        {
            // new group
            // test if it's one that we want
            //std::cout << "Group_" << XMLCH_2_STR(element->getAttribute(ATTR_gid)) << "_" << XMLCH_2_STR(element->getAttribute(ATTR_drseq)) << ".fa" << std::endl;
            char * c_group_name = tc(currentElement->getAttribute(ATTR_gid));
            std::string current_group_name = c_group_name;
            xr(&c_group_name);
            if (current_group_name == wantedGroup) 
            {
                // get the length of the direct repeat
                char * c_dr = tc(currentElement->getAttribute(ATTR_drseq));
                *directRepeat = c_dr;
                return currentElement;
            }
        }
        
    }
    
    // we should theoretically never get here but if the xml is bad then it might just happen
    // or if the user has put in a group that doesn't exist by mistake
    return NULL;
}

xercesc::DOMElement * crispr::XML::parseGroupForAssembly(xercesc::DOMElement* parentNode)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
       if( XMLString::equals(currentElement->getTagName(), TAG_assembly))
       {
           // assembly section
           // the child nodes will be the contigs
           return currentElement;
       }       
   }
    // if there is no assembly for this group
    return NULL;
} 

void crispr::XML::parseAssemblyForContigIds(xercesc::DOMElement* parentNode, std::set<std::string>& wantedContigs, std::list<std::string>& spacersForAssembly)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
       if( XMLString::equals(currentElement->getTagName(), TAG_contig))
       {
           // check to see if the current contig is one that we want
           char * c_current_contig = tc(currentElement->getAttribute(ATTR_cid));
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

void crispr::XML::getSpacerIdForAssembly(xercesc::DOMElement* parentNode, std::list<std::string>& spacersForAssembly)
{
    for (xercesc::DOMElement * currentElement = parentNode->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling())        
    {
        if( XMLString::equals(currentElement->getTagName(), TAG_cspacer))
        {
            char * c_cspacer = tc(currentElement->getAttribute(ATTR_spid));
            std::string str = c_cspacer;
            spacersForAssembly.push_back(str);
            xr(&c_cspacer);
        }
    }
}

DOMDocument * crispr::XML::setFileParser(const char * XMLFile)
{
    // Configure DOM parser.
    CX_FileParser->setValidationScheme( XercesDOMParser::Val_Never );
    CX_FileParser->setDoNamespaces( false );
    CX_FileParser->setDoSchema( false );
    CX_FileParser->setLoadExternalDTD( false );
    
    try
    {
        CX_FileParser->parse( XMLFile );
        return CX_FileParser->getDocument();        
    }
    catch( xercesc::XMLException& e ) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        XMLString::release( &message );
        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
    } catch (xercesc::DOMException& e) {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::stringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        XMLString::release( &message );
        throw crispr::xml_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(errBuf.str()).c_str());
    }
}

DOMElement * crispr::XML::createDOMDocument(std::string rootElement, std::string versionNumber, int& errorNumber )   
{
    XMLCh * core = tc("Core");
    DOMImplementation* impl =  DOMImplementationRegistry::getDOMImplementation(core);
    xr(&core);
    if (impl != NULL)
    {
        try
        {
            XMLCh * x_root_elem = tc(rootElement.c_str());
            CX_DocElem = impl->createDocument( 0, x_root_elem, 0);  
            xr(&x_root_elem);

            if (CX_DocElem != NULL) 
            {
                DOMElement* rootElem = CX_DocElem->getDocumentElement();
                XMLCh * x_version_num = tc(versionNumber.c_str());

                rootElem->setAttribute(ATTR_version, x_version_num);
                xr(&x_version_num);

                errorNumber = 0;
                return rootElem;
            }
        }
        catch (const OutOfMemoryException&)
        {
            XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
            errorNumber =  5;
        }
        catch (const DOMException& e)
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

DOMElement * crispr::XML::createDOMDocument(const char * rootElement, const char * versionNumber, int& errorNumber )   
{
    XMLCh * core = tc("Core");
    DOMImplementation* impl =  DOMImplementationRegistry::getDOMImplementation(core);
    xr(&core);
    if (impl != NULL)
    {
        try
        {
            XMLCh * x_root_elem = tc(rootElement);
            CX_DocElem = impl->createDocument(0, x_root_elem, 0);
            
            xr(&x_root_elem);
            
            if (CX_DocElem != NULL) 
            {
                DOMElement* rootElem = CX_DocElem->getDocumentElement();
                XMLCh * x_version_num = tc(versionNumber);
                
                rootElem->setAttribute(ATTR_version, x_version_num );
                
                xr(&x_version_num);
                errorNumber = 0;
                return rootElem;
            }
        }
        catch (const OutOfMemoryException&)
        {
            XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
            errorNumber =  5;
        }
        catch (const DOMException& e)
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

xercesc::DOMElement * crispr::XML::addMetaData(std::string notes, DOMElement * parentNode)
{
    DOMElement * meta_data_elem = CX_DocElem->createElement(TAG_metadata);
    DOMElement * notes_elem = CX_DocElem->createElement(TAG_notes);
    XMLCh * x_notes = tc(notes.c_str());
    DOMText * meta_data_notes = CX_DocElem->createTextNode(x_notes);
    xr(&x_notes);
    notes_elem->appendChild(meta_data_notes);
    meta_data_elem->appendChild(notes_elem);
    parentNode->appendChild(meta_data_elem);
    
    return meta_data_elem;
    
}
void crispr::XML::addFileToMetadata(std::string type, std::string url, DOMElement * parentNode)
{
    DOMElement * file = CX_DocElem->createElement(TAG_file);
    XMLCh * x_type = tc(type.c_str());
    XMLCh * x_url = tc(url.c_str());
    file->setAttribute(ATTR_type, x_type);
    file->setAttribute(ATTR_url, x_url);
    xr(&x_type);
    xr(&x_url);
    parentNode->appendChild(file);
}
xercesc::DOMElement * crispr::XML::addGroup(std::string& gID, std::string& drConsensus, DOMElement * parentNode)
{
    DOMElement * group = CX_DocElem->createElement(TAG_group);
    
    // Set the attributes of the group
    XMLCh * x_gID = tc(gID.c_str());
    XMLCh * x_drConsensus = tc(drConsensus.c_str());
    
    group->setAttribute(ATTR_gid, x_gID);
    group->setAttribute(ATTR_drseq, x_drConsensus);
    
    xr(&x_gID);
    xr(&x_drConsensus);
    // add the group to the parent (root element)
    parentNode->appendChild(group);
    return group;
}
xercesc::DOMElement * crispr::XML::addData(xercesc::DOMElement * parentNode)
{
    // create the data node with spacer and drs as child elements 
    DOMElement * data = CX_DocElem->createElement(TAG_data);
    DOMElement * drs = CX_DocElem->createElement(TAG_drs);
    DOMElement * spacers = CX_DocElem->createElement(TAG_spacers);
    data->appendChild(drs);
    data->appendChild(spacers);
    parentNode->appendChild(data);
    return data;
}
xercesc::DOMElement * crispr::XML::addAssembly(xercesc::DOMElement * parentNode)
{
    DOMElement * assembly = CX_DocElem->createElement(TAG_assembly);
    parentNode->appendChild(assembly);
    return assembly;
}
void crispr::XML::addDirectRepeat(std::string& drid, std::string& seq, DOMElement * parentNode)
{
    DOMElement * dr = CX_DocElem->createElement(TAG_dr);
    
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_drid = tc(drid.c_str());
    
    dr->setAttribute(ATTR_seq, x_seq);
    dr->setAttribute(ATTR_drid, x_drid);
    
    xr(&x_seq);
    xr(&x_drid);
    
    parentNode->appendChild(dr);
}
void crispr::XML::addSpacer(std::string& seq, std::string& spid, DOMElement * parentNode, std::string cov)
{
    DOMElement * sp = CX_DocElem->createElement(TAG_spacer);
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_cov = tc(cov.c_str());
    XMLCh * x_spid = tc(spid.c_str());
    
    sp->setAttribute(ATTR_seq, x_seq);
    sp->setAttribute(ATTR_spid, x_spid);
    sp->setAttribute(ATTR_cov, x_cov);
    
    xr(&x_seq);
    xr(&x_cov);
    xr(&x_spid);
    
    parentNode->appendChild(sp);
}
xercesc::DOMElement * crispr::XML::createFlankers(xercesc::DOMElement * parentNode)
{
    DOMElement * flankers = CX_DocElem->createElement(TAG_flankers);
    parentNode->appendChild(flankers);
    return flankers;
}
void crispr::XML::addFlanker(std::string& seq, std::string& flid, xercesc::DOMElement * parentNode)
{
    XMLCh * x_seq = tc(seq.c_str());
    XMLCh * x_flid = tc(flid.c_str());

    DOMElement * flanker = CX_DocElem->createElement(TAG_flanker);
    flanker->setAttribute(ATTR_seq, x_seq);
    flanker->setAttribute(ATTR_flid, x_flid);
    
    xr(&x_seq);
    xr(&x_flid);
    
    parentNode->appendChild(flanker);
}
xercesc::DOMElement * crispr::XML::addContig(std::string& cid, DOMElement * parentNode)
{
    DOMElement * contig = CX_DocElem->createElement(TAG_contig);
    
    XMLCh * x_cid = tc(cid.c_str());
    contig->setAttribute(ATTR_cid, x_cid);
    xr(&x_cid);
    parentNode->appendChild(contig);
    return contig;
}
void crispr::XML::createConsensus(std::string& concensus, xercesc::DOMElement * parentNode)
{
    DOMElement * concensus_elem = CX_DocElem->createElement(TAG_consensus);
    XMLCh * x_consensus = tc(concensus.c_str());
    DOMText * concensus_text = CX_DocElem->createTextNode(x_consensus);
    
    xr(&x_consensus);
    
    concensus_elem->appendChild(concensus_text);
    parentNode->appendChild(concensus_elem);
}
xercesc::DOMElement * crispr::XML::addSpacerToContig(std::string& spid, DOMElement * parentNode)
{
    DOMElement * cspacer = CX_DocElem->createElement(TAG_cspacer);
    XMLCh * x_spid = tc(spid.c_str());
    cspacer->setAttribute(ATTR_spid, x_spid);
    xr(&x_spid);
    parentNode->appendChild(cspacer);
    return cspacer;
}
xercesc::DOMElement * crispr::XML::createSpacers(std::string tag)
{
    XMLCh * x_tag = tc(tag.c_str());
    DOMElement * spacers = CX_DocElem->createElement(x_tag);
    xr(&x_tag);
    return spacers;
}

xercesc::DOMElement * crispr::XML::createFlankers(std::string tag)
{
    return createSpacers(tag);
}

void crispr::XML::addSpacer(std::string tag, std::string& spid, std::string& drid, std::string& drconf, DOMElement * parentNode)
{
    XMLCh * x_tag = tc(tag.c_str());
    XMLCh * x_drid = tc(drid.c_str());
    XMLCh * x_drconf = tc(drconf.c_str());
    XMLCh * x_spid = tc(spid.c_str());

    DOMElement * fs = CX_DocElem->createElement(x_tag);
    fs->setAttribute(ATTR_drid, x_drid);
    fs->setAttribute(ATTR_drconf, x_drconf);
    fs->setAttribute(ATTR_spid, x_spid);
    
    xr(&x_tag);
    xr(&x_drid);
    xr(&x_drconf);
    xr(&x_spid);
    
    parentNode->appendChild(fs);
}
void crispr::XML::addFlanker(std::string tag, std::string& flid, std::string& drconf, std::string& directjoin, xercesc::DOMElement * parentNode)
{
    XMLCh * x_tag = tc(tag.c_str());
    XMLCh * x_flid = tc(flid.c_str());
    XMLCh * x_drconf = tc(drconf.c_str());
    XMLCh * x_directjoin = tc(directjoin.c_str());
    
    DOMElement * bf = CX_DocElem->createElement(x_tag);
    bf->setAttribute(ATTR_flid, x_flid);
    bf->setAttribute(ATTR_drconf, x_drconf);
    bf->setAttribute(ATTR_directjoin, x_directjoin);
    
    xr(&x_tag);
    xr(&x_flid);
    xr(&x_drconf);
    xr(&x_directjoin);
    
    parentNode->appendChild(bf);
}

//
 // File IO / Printing
//

bool crispr::XML::printDOMToFile(std::string outFileName )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {chLatin_L, chLatin_S, chNull};
        DOMImplementation *impl          = DOMImplementationRegistry::getDOMImplementation(tempStr);
        DOMLSSerializer   *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();
        DOMLSOutput       *theOutputDesc = ((DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);

        DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        XMLFormatTarget *myFormTarget;
        myFormTarget = new LocalFileFormatTarget(outFileName.c_str());
        //myFormTarget=new StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);

        theSerializer->write(CX_DocElem, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (XMLException& e)
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

bool crispr::XML::printDOMToScreen(void )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {chLatin_L, chLatin_S, chNull};
        DOMImplementation *impl          = DOMImplementationRegistry::getDOMImplementation(tempStr);
        DOMLSSerializer   *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();
        DOMLSOutput       *theOutputDesc = ((DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        XMLFormatTarget *myFormTarget;
        myFormTarget=new StdOutFormatTarget();
        
        theOutputDesc->setByteStream(myFormTarget);
        
        theSerializer->write(CX_DocElem, theOutputDesc);
        
        theOutputDesc->release();
        theSerializer->release();
        
        //
        // Filter, formatTarget and error handler
        // are NOT owned by the serializer.
        //
        delete myFormTarget;
        retval = true;
        
    }
    catch (const OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (XMLException& e)
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

bool crispr::XML::printDOMToFile(std::string outFileName, DOMDocument * docDOM )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {chLatin_L, chLatin_S, chNull};
        DOMImplementation *impl          = DOMImplementationRegistry::getDOMImplementation(tempStr);
        DOMLSSerializer   *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();
        DOMLSOutput       *theOutputDesc = ((DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");
        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        
        DOMConfiguration* serializerConfig = theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        XMLFormatTarget *myFormTarget;
        myFormTarget = new LocalFileFormatTarget(outFileName.c_str());
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
    catch (const OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (XMLException& e)
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

bool crispr::XML::printDOMToScreen(DOMDocument * domDoc )
{
    bool retval;
    
    try
    {
        // get a serializer, an instance of DOMLSSerializer
        XMLCh tempStr[3] = {chLatin_L, chLatin_S, chNull};
        DOMImplementation *impl          = DOMImplementationRegistry::getDOMImplementation(tempStr);
        DOMLSSerializer   *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer();
        DOMLSOutput       *theOutputDesc = ((DOMImplementationLS*)impl)->createLSOutput();
        
        // set user specified output encoding
        XMLCh * x_encoding = tc("ISO8859-1");

        theOutputDesc->setEncoding(x_encoding);
        xr(&x_encoding);
        
        
        DOMConfiguration* serializerConfig=theSerializer->getDomConfig();
        
        // set feature if the serializer supports the feature/mode
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTSplitCdataSections, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTSplitCdataSections, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
            serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
        
        if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTBOM, false))
            serializerConfig->setParameter(XMLUni::fgDOMWRTBOM, false);
        
        //
        // Plug in a format target to receive the resultant
        // XML stream from the serializer.
        //
        // StdOutFormatTarget prints the resultant XML stream
        // to stdout once it receives any thing from the serializer.
        //
        XMLFormatTarget *myFormTarget;
        myFormTarget=new StdOutFormatTarget();
        
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
    catch (const OutOfMemoryException&)
    {
        XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
        retval = false;
    }
    catch (XMLException& e)
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


