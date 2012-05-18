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
#include <set>
#include "kseq.h"
#include <fstream>

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


} assemblyOptions;

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
