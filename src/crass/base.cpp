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

#include "base.h"
using namespace crispr::xml;


base::base()
{
    try  
    {
        // Initialize Xerces infrastructure   
        init();
        // transcode all the member variables
        alloc();
    }
    catch( xercesc::XMLException& e )
    {
        char * message = xercesc::XMLString::transcode( e.getMessage() );
        std::cerr << "XML toolkit initialization error: " << message << std::endl;
        xercesc::XMLString::release( &message );
        // throw exception here to return ERROR_XERCES_INIT
    }
}


base::~base(void)
{
    // Free memory
    dealloc();
    // Terminate Xerces
    release();
}
void base::init(void) {
    xercesc::XMLPlatformUtils::Initialize();
}
void base::alloc(void) {
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

void base::dealloc(void) {
    // grep ELEMENT crass.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLString::release( \&TAG_\1 );%" | sort | uniq
    // grep ATTLIST crass.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLString::release( \&ATTR_\1 );%" | sort | uniq
    try {
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
    } catch(...) {
        throw crispr::xml_exception (__FILE__,
                                     __LINE__,
                                     __PRETTY_FUNCTION__,
                                     "Unknown error occurred");
    }
    
}

void base::release (void) {
    xercesc::XMLPlatformUtils::Terminate();  // Terminate after release of memory
}

