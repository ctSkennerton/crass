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

#ifndef BASE_H
#define BASE_H

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
#include <map>

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif
#include "Exception.h"


#define tc(buf) xercesc::XMLString::transcode(buf)
// same as above but more verbose
#define transcode_xmlch(buf) xercesc::XMLString::transcode(buf)
#define xr(buf) xercesc::XMLString::release(buf)
// same as above but more verbose
#define free_xmlch(buf) xercesc::XMLString::release(buf)
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
            void release(void);
            //
            // Generic get
            //
            // grep ATTLIST crispr-1.1.dtd | perl -ne 's/[^ ]* [^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * attr_$_(void) { return ATTR_$original; }\n";' | sort | uniq
            // grep ELEMENT crispr-1.1.dtd | perl -ne 's/[^ ]* ([^ ]*) .*/\1/ ;chomp; my $original = $_; s/\b(\w)/\U$1/g; print "inline XMLCh * tag_$_(void) { return TAG_$original; }\n";' |sort | uniq
            
            
            /** Accessor to the 'accession' attribute   
             *  @return XMLCh * of 'accession'  
             */
            inline XMLCh * attr_Accession(void) { return ATTR_accession; }
            
            /** Accessor to the 'cid' attribute   
             *  @return XMLCh * of 'cid'  
             */
            inline XMLCh * attr_Cid(void) { return ATTR_cid; }
            
            /** Accessor to the 'confcnt' attribute   
             *  @return XMLCh * of 'confcnt'  
             */
            inline XMLCh * attr_Confcnt(void) { return ATTR_confcnt; }
            
            /** Accessor to the 'cov' attribute   
             *  @return XMLCh * of 'cov'  
             */
            inline XMLCh * attr_Cov(void) { return ATTR_cov; }
            
            /** Accessor to the 'directjoin' attribute   
             *  @return XMLCh * of 'directjoin'  
             */
            inline XMLCh * attr_Directjoin(void) { return ATTR_directjoin; }
            
            /** Accessor to the drconf attribute
             *   @return  XMLCh * of drconf
             */
            inline XMLCh * attr_Drconf(void) { return ATTR_drconf; }
            
            /** Accessor to the drid attribute
             *  @return  XMLCh * of drid
             */
            inline XMLCh * attr_Drid(void) { return ATTR_drid; }
            
            /** Accessor to the drseq attribute
             *  @return  XMLCh * of drseq
             */
            inline XMLCh * attr_Drseq(void) { return ATTR_drseq; }
            
            /** Accessor to the flid attribute
             *  @return  XMLCh * of flid
             */
            inline XMLCh * attr_Flid(void) { return ATTR_flid; }
            
            /** Accessor to the gid attribute
             *  @return  XMLCh * of gid
             */
            inline XMLCh * attr_Gid(void) { return ATTR_gid; }
            
            /** Accessor to the seq attribute
             *  @return  XMLCh * of seq
             */
            inline XMLCh * attr_Seq(void) { return ATTR_seq; }
            
            /** Accessor to the soid attribute
             *  @return  XMLCh * of soid
             */
            inline XMLCh * attr_Soid(void) { return ATTR_soid; }
            
            /** Accessor to the spid attribute
             *  @return  XMLCh * of spid
             */
            inline XMLCh * attr_Spid(void) { return ATTR_spid; }
            
            /** Accessor to the totcnt attribute
             *  @return  XMLCh * of totcnt
             */
            inline XMLCh * attr_Totcnt(void) { return ATTR_totcnt; }
            
            /** Accessor to the type attribute
             *  @return  XMLCh * of type
             */
            inline XMLCh * attr_Type(void) { return ATTR_type; }
            
            /** Accessor to the url attribute
             *  @return  XMLCh * of url
             */
            inline XMLCh * attr_Url(void) { return ATTR_url; }
            
            /** Accessor to the version attribute
             *  @return  XMLCh * of version
             */
            inline XMLCh * attr_Version(void) { return ATTR_version; }
            
            /** Accessor to the assembly element
             *  @return  XMLCh * of assembly
             */
            inline XMLCh * tag_Assembly(void) { return TAG_assembly; }
            
            /** Accessor to the bf element
             *  @return  XMLCh * of bf
             */
            inline XMLCh * tag_Bf(void) { return TAG_bf; }
            
            /** Accessor to the bflankers element
             *  @return  XMLCh * of bflankers
             */
            inline XMLCh * tag_Bflankers(void) { return TAG_bflankers; }
            
            /** Accessor to the bs element
             *  @return  XMLCh * of bs
             */
            inline XMLCh * tag_Bs(void) { return TAG_bs; }
            
            /** Accessor to the bspacers element
             *  @return  XMLCh * of bspacers
             */
            inline XMLCh * tag_Bspacers(void) { return TAG_bspacers; }
            
            /** Accessor to the command element
             *  @return  XMLCh * of command
             */
            inline XMLCh * tag_Command(void) { return TAG_command; }
            
            /** Accessor to the consensus element
             *  @return  XMLCh * of consensus
             */
            inline XMLCh * tag_Consensus(void) { return TAG_consensus; }
            
            /** Accessor to the contig element
             *  @return  XMLCh * of contig
             */
            inline XMLCh * tag_Contig(void) { return TAG_contig; }
            
            /** Accessor to the crispr element
             *  @return  XMLCh * of crispr
             */
            inline XMLCh * tag_Crispr(void) { return TAG_crispr; }
            
            /** Accessor to the cspacer element
             *  @return  XMLCh * of cspacer
             */
            inline XMLCh * tag_Cspacer(void) { return TAG_cspacer; }
            
            /** Accessor to the data element
             *  @return  XMLCh * of data
             */
            inline XMLCh * tag_Data(void) { return TAG_data; }
            
            /** Accessor to the dr element
             *  @return  XMLCh * of dr
             */
            inline XMLCh * tag_Dr(void) { return TAG_dr; }
            
            /** Accessor to the drs element
             *  @return  XMLCh * of drs
             */
            inline XMLCh * tag_Drs(void) { return TAG_drs; }
            
            /** Accessor to the epos element
             *  @return  XMLCh * of epos
             */
            inline XMLCh * tag_Epos(void) { return TAG_epos; }
            
            /** Accessor to the ff element
             *  @return  XMLCh * of ff
             */
            inline XMLCh * tag_Ff(void) { return TAG_ff; }
            
            /** Accessor to the fflankers element
             *  @return  XMLCh * of fflankers
             */
            inline XMLCh * tag_Fflankers(void) { return TAG_fflankers; }
            
            /** Accessor to the file element
             *  @return  XMLCh * of file
             */
            inline XMLCh * tag_File(void) { return TAG_file; }
            
            /** Accessor to the flanker element
             *  @return  XMLCh * of flanker
             */
            inline XMLCh * tag_Flanker(void) { return TAG_flanker; }
            
            /** Accessor to the flankers element
             *  @return  XMLCh * of flankers
             */
            inline XMLCh * tag_Flankers(void) { return TAG_flankers; }
            
            /** Accessor to the fs element
             *  @return  XMLCh * of fs
             */
            inline XMLCh * tag_Fs(void) { return TAG_fs; }
            
            /** Accessor to the fspacers element
             *  @return  XMLCh * of fspacers
             */
            inline XMLCh * tag_Fspacers(void) { return TAG_fspacers; }
            
            /** Accessor to the group element
             *  @return  XMLCh * of group
             */
            inline XMLCh * tag_Group(void) { return TAG_group; }
            
            /** Accessor to the metadata element
             *  @return  XMLCh * of metadata
             */
            inline XMLCh * tag_Metadata(void) { return TAG_metadata; }
            
            /** Accessor to the name element
             *  @return  XMLCh * of name
             */
            inline XMLCh * tag_Name(void) { return TAG_name; }
            
            /** Accessor to the notes element
             *  @return  XMLCh * of notes
             */
            inline XMLCh * tag_Notes(void) { return TAG_notes; }
            
            /** Accessor to the program element
             *  @return  XMLCh * of program
             */
            inline XMLCh * tag_Program(void) { return TAG_program; }
            
            /** Accessor to the source element
             *  @return  XMLCh * of source
             */
            inline XMLCh * tag_Source(void) { return TAG_source; }
            
            /** Accessor to the sources element
             *  @return  XMLCh * of sources
             */
            inline XMLCh * tag_Sources(void) { return TAG_sources; }
            
            /** Accessor to the spos element
             *  @return  XMLCh * of spos
             */
            inline XMLCh * tag_Spos(void) { return TAG_spos; }
            
            /** Accessor to the spacer element
             *  @return  XMLCh * of spacer
             */
            inline XMLCh * tag_Spacer(void) { return TAG_spacer; }
            
            /** Accessor to the spacers element
             *  @return  XMLCh * of spacers
             */
            inline XMLCh * tag_Spacers(void) { return TAG_spacers; }
            
            /** Accessor to the version element
             *  @return  XMLCh * of version
             */
            inline XMLCh * tag_Version(void) { return TAG_version; }
            
            
            
        private:            
            // grep ATTLIST crispr-1.1.dtd | sed -e "s%[^ ]* [^ ]* \([^ ]*\) .*%XMLCh\* ATTR_\1;%" | sort | uniq
            // grep ELEMENT crass-1.1.dtd | sed -e "s%[^ ]* \([^ ]*\) .*%XMLCh\* TAG_\1;%" | sort | uniq
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
    }
}

#endif 