/*
 *  AssemblyWrapper.h is part of the crass project
 *  
 *  Created by Connor Skennerton on 27/10/11.
 *  Copyright 2011, 2012 Connor Skennerton and Michael Imelfort. All rights reserved. 
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

#ifndef crass_AssemblyWrapper_h
#define crass_AssemblyWrapper_h

#include <string>
#include <vector>
#include <map>
#include <set>
#include "kseq.h"
#include <fstream>
#include "reader.h"

enum ASSEMBLERS {
    unset,
    velvet,
    cap3
    };

typedef struct{
    std::string xmlFileName;            // name of the crass xml file
    int         group;                  // group id to parse
    int         logLevel;               // log level
    bool        pairedEnd;              // does the user want a paired end assembly
    std::string inputDirName;           // location of the input directory
    std::string outputDirName;          // location of the output directory
    std::string segments;               // a comma separated list of segments to assemble
    bool        logToScreen;            // does the user want logging info printed to screen
    int         insertSize;             // the insert size for the paired end assembly
    ASSEMBLERS  assembler;              // the assembler that the user wants

} assemblyOptions;

typedef std::vector<std::string> IDVector;
typedef std::map<std::string, IDVector > Spacer2SourceMap;
typedef std::map<std::string, std::string > XMLIDMap;
typedef std::set<std::string> StringSet;

class CrisprParser : public crispr::xml::reader {
    
    
public:
    //
    // Working functions
    //
    void parseXMLFile(std::string XMLFile, 
                      std::string& wantedGroup, 
                      std::string& directRepeat,
                      StringSet& wantedContigs,
                      StringSet& wantedSpacers
                      );
    
    xercesc::DOMElement * getWantedGroupFromRoot(xercesc::DOMElement * currentElement, 
                                                 std::string& wantedGroup, 
                                                 std::string&  directRepeat
                                                 );
    
    
    xercesc::DOMElement * parseGroupForAssembly(xercesc::DOMElement* currentElement
                                                );
    
    void parseAssemblyForContigIds(xercesc::DOMElement* parentNode, 
                                   StringSet& wantedReads, 
                                   Spacer2SourceMap& spacersForAssembly,
                                   XMLIDMap& source2acc,
                                   StringSet& wantedContigs
                                   );
    /** Iterate through all spacers for a contig and collecting source accessions  
     *  @param parentNode The contig node from the DOM tree 
     *  @param wantedReads a container to store the source accessions in
     *  @param spacersForAssembly A container containing a list of all sources for each spacer
     *  @param source2acc A container containing a mapping between source IDs and source accessions 
     */
    void getSourceIdForAssembly(xercesc::DOMElement* parentNode, 
                                StringSet& wantedReads,
                                Spacer2SourceMap& spacersForAssembly,
                                XMLIDMap& source2acc
                                );
    
    /** Get a list of sources for a group 
     *  @param container An associative container to write the sources to. 
     *  The container must have both key and value types as XMLCh *
     *  The container must overload reference operator[]
     *  @param parentNode The parent node for the sources for iteration
     */
    template <class C>
    void getSourcesForGroup(C& container, xercesc::DOMElement * parentNode) {
        
        if (xercesc::XMLString::equals(tag_Sources(), parentNode->getTagName())) {
            for (xercesc::DOMElement * currentNode = parentNode->getFirstElementChild(); 
                 currentNode != NULL; 
                 currentNode = currentNode->getNextElementSibling()) {
                
                
                char * current_source_id = tc(currentNode->getAttribute(attr_Soid()));
                char * current_source_acc = tc(currentNode->getAttribute(attr_Accession()));
                container[current_source_id] = current_source_acc;
                xr(&current_source_id);
                xr(&current_source_acc);
            }
        } else {
            throw crispr::xml_exception(__FILE__,
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        "Parent not 'sources' tag");
        }
        
        
    }
    /** Get a list of source identifiers for all spacers
     *  @param container An associative container that maps a single
     *  key (the spacer) to a list of sources
     *  @param parentNode The parent node for the spacer tags
     */
    template <class C>
    void mapSacersToSourceID(C& container, xercesc::DOMElement * parentNode) {
        // each spacer
        for (xercesc::DOMElement * currentNode = parentNode->getFirstElementChild(); 
             currentNode != NULL; 
             currentNode = currentNode->getNextElementSibling()) {
            
            char * spid = tc(currentNode->getAttribute(attr_Spid()));
            // each source
            for (xercesc::DOMElement * sp_source = currentNode->getFirstElementChild(); 
                 sp_source != NULL; 
                 sp_source = sp_source->getNextElementSibling()) {
                
                char * soid = tc(sp_source->getAttribute(attr_Soid()));
                container[spid].push_back(soid);
                xr(&soid);
            }
            xr(&spid);
        }
    }
};


void assemblyVersionInfo(void); 

void assemblyUsage(void);

int processAssemblyOptions(int argc, char * argv[], assemblyOptions& opts);

inline int calculateOverlapLength(int drLength){return drLength + 8;}

void parseSegmentString(std::string& segmentString, std::set<std::string>& segments);

void generateTmpAssemblyFile(std::string fileName, std::set<std::string>& wantedContigs, assemblyOptions& opts, std::string& tmpFileName);

int velvetWrapper(int hashLength, assemblyOptions& opts, std::string& tmpFileName);

int capWrapper(int overlapLength, assemblyOptions& opts, std::string& tmpFileName);

int assemblyMain(int argc, char * argv[]);


#endif
