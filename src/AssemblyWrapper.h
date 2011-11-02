/*
 *  AssemblyWrapper.h is part of the crass project
 *  
 *  Created by Connor Skennerton on 27/10/11.
 *  Copyright 2011 Connor Skennerton and Michael Imelfort. All rights reserved. 
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

void assemblyUsage(void);

int processAssemblyOptions(int argc, char * argv[], assemblyOptions& opts);

int calculateOverlapLength(int group, std::string& inputDir);

void parseSegmentString(std::string& segmentString);

void velvetWrapper(int hashLength, std::string& inputDir);

void capWrapper(std::string& inputFile, std::string& inputDir);

void assemblyMain(int argc, char * argv[]);


// templated function to split a string on delimeters and return it in a container of class T
template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", const bool trimEmpty = false)
{
    std::string::size_type pos, lastPos = 0;
    while(true)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos)
        {
            pos = str.length();
            
            if(pos != lastPos || !trimEmpty)
                tokens.push_back( typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos - lastPos ));
            
            break;
        }
        else
        {
            if(pos != lastPos || !trimEmpty)
                tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos ));
        }
        
        lastPos = pos + 1;
    }
}

#endif
