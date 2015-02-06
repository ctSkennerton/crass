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

#ifndef WRITER_H
#define WRITER_H 
#include "base.h"

namespace crispr {
    namespace xml {
        class writer : virtual public base {
            
            //members
            xercesc::DOMDocument * XW_DocElem;
            int XW_CurrentSourceId;
            
        public:
            
            //constructor/destructor
            writer();
            ~writer();
            
            /** Create an in-memory representation of a crispr file   
             *  @param rootElement Name for the root element 
             *  @param versionNumber version for the crispr file to have
             *  @param errorNumber an integer for saving any error codes produced during document creation
             *  @return The xercesc::DOMElement for the root node or NULL on failure
             */
            xercesc::DOMElement * createDOMDocument(std::string rootElement, std::string versionNumber, int& errorNumber );
            
            xercesc::DOMElement * createDOMDocument(const char * rootElement, const char * versionNumber, int& errorNumber );
            
            /** add a 'metadata' tag to 'group'   
             *  @param parentNode the xercesc::DOMElement of the 'group' tag 
             *  @return  the xercesc::DOMElement of the 'metadata' tag
             */
            xercesc::DOMElement * addMetaData(xercesc::DOMElement * parentNode);
            
            /** add a 'metadata' tag to 'group' with notes  
             *  @param parentNode the xercesc::DOMElement of the 'group' tag
             *  @param notes freeform notes that the user can add in 
             *  @return  the xercesc::DOMElement of the 'metadata' tag
             */
            xercesc::DOMElement * addMetaData(std::string notes, xercesc::DOMElement * parentNode);
            
            /** add in files referenced by this group  
             *  @param type The type of reference file.  Must be one of image|sequence|log|data
             *  @param url path to the file
             *  @param parentNode xercesc::DOMElement of the 'metadata' tag
             */
            void addFileToMetadata(std::string type, std::string url, xercesc::DOMElement * parentNode);
            
            /** add extra notes to the 'metadata'  
             *  @param notes freeform notes added by the user
             *  @param parentNode xercesc::DOMElement of the 'metadata' tag
             */
            void addNotesToMetadata(std::string notes, xercesc::DOMElement * parentNode);
            
            /** add a 'group' to the root element  
             *  @param gID The unique group identifier
             *  @param drConsensus The direct repeat concensus sequence for this group
             *  @param parentNode xercesc::DOMElement of the root element
             */
            xercesc::DOMElement * addGroup(std::string& gID, std::string& drConsensus, xercesc::DOMElement * parentNode);
            
            /** add the 'data' to 'group'. This automatically creates the child nodes for the 'sources' 'drs' and 'spacers' tags  
             *  @param parentNode the 'group' xercesc::DOMElement
             *  @return The xercesc::DOMElement of the 'data' tag   
             */
            xercesc::DOMElement * addData(xercesc::DOMElement * parentNode);
            
            /** add the 'assembly' to 'group'   
             *  @param parentNode the 'group' xercesc::DOMElement
             *  @return The xercesc::DOMElement of the 'assembly' tag   
             */
            xercesc::DOMElement * addAssembly(xercesc::DOMElement * parentNode);
            
            /** add a direct repeat 'dr' tag to the 'drs' tag   
             *  @param drid The unique DR id for this group
             *  @param seq The sequence of the direct repeat in its lowest lexicographical form
             *  @param parentNode The xercesc::DOMElement of the 'drs' tag 
             */
            void addDirectRepeat(std::string& drid, std::string& seq, xercesc::DOMElement * parentNode);
            
            /** add a 'spcaer' tag to the 'spacers'   
             *  @param seq sequence of the spacer in the form that would appear if the direct repeat is in its lowest lexicographical form
             *  @param spid  the unique spacer identifier for this group
             *  @param parentNode the xercesc::DOMElement of the 'spacers' tag
             *  @param cov The coverage of the spacer. Defaults to zero
             *  @return the xercesc::DOMElement for the 'spacer' tag
             */
            xercesc::DOMElement * addSpacer(std::string& seq, std::string& spid, xercesc::DOMElement * parentNode, std::string cov = "0" );
            
            /** create a 'flankers' tag in 'data'   
             *  @param parentNode the xercesc::DOMElement of the 'data' tag
             *  @return the xercesc::DOMElement of the 'flankers' tag
             */
            xercesc::DOMElement * createFlankers(xercesc::DOMElement * parentNode);
            
            /** add a 'flanker' tag to the 'flankers' tag   
             *  @param seq sequence of the flanker in the form that would appear if the direct repeat is in its lowest lexicographical form
             *  @param flid the unique flanker identifier for this group
             *  @param parentNode the xercesc::DOMElement of the 'flankers' tag 
             *  @return the xercesc::DOMElement of the 'flanker' tag   
             */
            xercesc::DOMElement * addFlanker(std::string& seq, std::string& flid, xercesc::DOMElement * parentNode);
            
            /** add a 'contig' to an 'assembly'   
             *  @param cid a unique contig id for this group
             *  @param parentNode the xercesc::DOMElement of the 'assembly' tag
             *  @return the xercesc::DOMElement of the 'contig' tag   
             */
            xercesc::DOMElement * addContig(std::string& cid, xercesc::DOMElement * parentNode);
            
            /** add the concensus sequence to the contig   
             *  @param concensus A string of DNA characters representing the joined sequence of all the spacers and direct repeats for this contig
             *  @param parentNode the xercesc::DOMElement of the 'contig' tag  
             */
            void createConsensus(std::string& concensus, xercesc::DOMElement * parentNode);
            
            /** Add a 'cspacer' to a 'contig'   
             *  @param spid Unique identifier for the spacer.  should be the same as listed in the 'spacers' tag in 'data'
             *  @param parentNode the xercesc::DOMElement of the 'contig' tag
             *  @return the xercesc::DOMElement of the 'cspacer' tag  
             */
            xercesc::DOMElement * addSpacerToContig(std::string& spid, xercesc::DOMElement * parentNode);
            
            /** use to create a backward spacers ('bspacers') or forward spacers ('fspacers') child elements of 'cspacer'   
             *  @param tag A string of either "bspacers" or "fspacers"
             *  @return the xercesc::DOMElement of either "bspacers" or "fspacers"  
             */
            xercesc::DOMElement * createSpacers(std::string tag);
            
            /** A convienience method that calls createSpacers, Used only for making the code more explicit   
             *  @param tag A string of either "bflankers" or "fflankers"
             *  @return the xercesc::DOMElement of either "bflankers" or "fflankers"  
             */
            xercesc::DOMElement * createFlankers(std::string tag);
            
            /** add either a backward spacer (bs) or forward spacer (fs) tags to a cspacer tag   
             *  @param tag A string that describes the type of association ('bs' or 'fs')
             *  @param spid The unique spacer id for the link.  Should be the same as an spid listed in 'spacers' tag inside 'data'
             *  @param drid The unique direct repeat id for the DR that sits inbetween the two spacers
             *  @param drconf
             *  @param parentNode The xercesc::DOMElement corresponding to one of bspacers, fspacers  
             */
            void addSpacer(std::string tag, std::string& spid, std::string& drid, std::string& drconf, xercesc::DOMElement * parentNode);
            
            /** add either a backward flanker (bf) or forward flanker (ff) tags to cspacer   
             *  @param tag A string that describes the type of association ('bf' or 'ff')
             *  @param flid The unique flanker id for the link.  Should be the same as an spid listed in 'fankers' tag inside 'data'
             *  @param drconf 
             *  @param directjoin
             *  @param parentNode The xercesc::DOMElement corresponding to one of bflankers, fflankers  
             */
            void addFlanker(std::string tag, std::string& flid, std::string& drconf, std::string& directjoin, xercesc::DOMElement * parentNode);
            
            /** create the sources tag for a group    
             *  @param parentNode The xercesc::DOMElement for the 'group' tag
             *  @return The xercesc::DOMElement for the 'sources' tag  
             */
            xercesc::DOMElement * addSources(xercesc::DOMElement * parentNode);
            
            /** add a 'source' tag for the sources   
             *  @param accession The accession for the source as it would appear in a fasta file
             *  @param soid The unique source identifier for this source.  Ideally this should be short, containing sharacters and
             *              numbers for future reference within the cirpsr file.
             *  @param parentNode The xercesc::DOMElement of the 'sources' tag
             *  @return The xercesc::DOMElement of teh 'source' tag  
             */
            xercesc::DOMElement * addSource(std::string accession, std::string soid, xercesc::DOMElement * parentNode);
            
            /** add a source tag for a spacer   
             *  @param soid The unique source identifier for this source. Sould be the same as listed in the 'sources' inside 'data'
             *  @param parentNode The xercesc::DOMElement of the 'spacer' that was fould inside this source
             *  @return The xercesc::DOMElement of the 'source' tag  
             */
            xercesc::DOMElement * addSpacerSource(std::string soid, xercesc::DOMElement * parentNode);
            
            /** add start and end positions for 'source' in 'spacer'   
             *  @param start The start position in the source for the spacer
             *  @param end The end position in the source for the spacer
             *  @param parentNode The 'source' tag to the spacer  
             */
            void addStartAndEndPos(std::string start, std::string end, xercesc::DOMElement * parentNode);
            
            /** add a 'program' tag to 'metadata'   
             *  @param parentNode The xercesc::DOMElement of the 'metadata' tag
             *  @return The xercesc::DOMElement of the 'program' tag  
             */
            xercesc::DOMElement * addProgram(xercesc::DOMElement * parentNode);
            
            /** add a 'name' tag to 'program'   
             *  @param progName The name of the program that called this CRISPR
             *  @param parentNode The xercesc::DOMElement of the 'program' tag
             */
            void addProgName(std::string progName, xercesc::DOMElement * parentNode);
            
            /** add a 'version' tag to 'program'   
             *  @param progVersion The version of the program that called this CRISPR
             *  @param parentNode The xercesc::DOMElement of the 'program' tag
             */
            void addProgVersion(std::string progVersion, xercesc::DOMElement * parentNode);
            
            /** add a 'command' tag to 'program'   
             *  @param progCommand The command line options used to call this CRISPR
             *  @param parentNode The xercesc::DOMElement of the 'program' tag 
             */
            void addProgCommand(std::string progCommand, xercesc::DOMElement * parentNode);
            
            /** print the current document to file   
             *  @param outFileName The name of the output file
             */
            bool printDOMToFile(std::string outFileName );
            
            /** print the current document to screen   
             */
            bool printDOMToScreen(void);
            
            bool printDOMToFile(std::string outFileName, xercesc::DOMDocument * domDoc );
            
            bool printDOMToScreen( xercesc::DOMDocument * domDoc);
            
            /** convienience method to return the root element of the current document   
             *  @return The xercesc::DOMElement for the root ('crispr') tag  
             */
            inline xercesc::DOMElement * getRootElement(void)
            {
                return XW_DocElem->getDocumentElement();
            }
            
            /** convienience method to return the current document   
             *  @return The xercesc::DOMDocument for this writer  
             */
            inline xercesc::DOMDocument * getDocumentObj(void)
            {
                return XW_DocElem;
            }
            
            
        };
    }
}

#endif