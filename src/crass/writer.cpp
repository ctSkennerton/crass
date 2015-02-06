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

#include "writer.h"


crispr::xml::writer::writer() {
    XW_DocElem = NULL;
}

crispr::xml::writer::~writer() {
    delete XW_DocElem;
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
    xercesc::DOMElement * data = XW_DocElem->createElement(tag_Data());
    xercesc::DOMElement * sources = XW_DocElem->createElement(tag_Sources());
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
