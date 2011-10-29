/*
 *  AssemblyWrapper.cpp is part of the crass project
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

// System Includes
#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <cstring>

// Local Includes
#include "AssemblyWrapper.h"
#include "config.h"
#include "crassDefines.h"
#include "StlExt.h"

void assemblyUsage(void)
{
    std::cout<<"Usage: "PACKAGE_NAME<<" assemble ASSEMBLER -g INT -s LIST -i INDIR [options]"<<std::endl;
    std::cout<<"\twhere ASSEMBLER is one of the assembly algorithms listed below:"<<std::endl<<std::endl;
#ifdef HAVE_VELVET
    std::cout<<"\tvelvet"<<std::endl;
#endif
    
#ifdef HAVE_CAP3
    std::cout<<"\tcap3"<<std::endl;
#endif
    std::cout<<std::endl;
    std::cout<< "-h --help                    This help message"<<std::endl;
    std::cout<< "-l --logLevel        <INT>   Output a log file and set a log level [1 - "<<CRASS_DEF_MAX_LOGGING<<"]"<<std::endl;
    std::cout<< "-V --version                 Program and version information"<<std::endl;
    std::cout<< "--logToScreen                Print the logging information to screen rather than a file"<<std::endl;
    std::cout<< "-g --group           <INT>   ID of the group that you want to assemble"<<std::endl;
    std::cout<< "-s --segments        <LIST>  A comma separated list of numbered segments to assemble from the specified group"<<std::endl;
    std::cout<< "-i --indir           <DIR>   input directory for the assembly [default: .]"<<std::endl;


}

int calculateOverlapLength(int group, std::string& inputDir)
{
    int olap = 1;
    //go to the input directory and look for the group required
    
    // return the length of the DR plus X number of bases
    return olap;
}


void parseSegmentString(std::string& segmentString)
{
    // split the segment id string into the individual pieces
    std::list<std::string> segments;
    tokenize(segmentString, segments, ",");
    
    // create a output stream for the tmp file
    std::ofstream combined_file;
    std::stringstream ss;
    ss<<PACKAGE_NAME<<"_tmp.fa";
    combined_file.open((ss.str()).c_str());
    
    if (combined_file.good()) {
        
        // iterate through the segments
        std::list<std::string>::iterator seg_iter = segments.begin();
        while (seg_iter != segments.end()) {
            
            // check for the file that should corespond to each token
            std::string curr_seg = "";
            
            // open the segment file
            std::ifstream curr_segment_file;
            curr_segment_file.open(curr_seg.c_str());
            
            if (curr_segment_file.good()) {
                
                // cat the files together
                char c;
                while (curr_segment_file >> c) combined_file << c;
                
                curr_segment_file.close();
                
            } else {
                // print error
            }
            seg_iter++;
        }
    } else {
        //TODO: print error
    }
}


void velvetWrapper( int hashLength, std::string& inputDir)
{
    // get the hash length from the DR
    
    
    // create the command string for velvet
    std::string h_cmd = "velveth " + inputDir +" " + to_string(hashLength) + " " + inputDir + "tmp.fa";
    std::string g_cmd = "velvetg " + inputDir ; 
    // execute velvet
    
    int h_exit = system(h_cmd.c_str());
    
    if (h_exit) {
        // print error
    }
    int g_exit = system(g_cmd.c_str());
    if (g_exit) {
        // print error
    }
    
}

void capWrapper(std::string& inputFile, std::string& inputDir)
{
    
}


void assemblyMain(int argc, char * argv[])
{
    if (argc == 1) {
        // print assembly wrapper help
        assemblyUsage();
        exit(1);
    }
    
#ifdef HAVE_VELVET
    if (strcmp(argv[1], "velvet") == 0) {
        // the user wants velvet
        
        // process options
        
        //pass everything into the velvet wrapper
    } 
#endif
    
#ifdef HAVE_CAP3
    if(strcmp(argv[1], "cap3") == 0){
        // the user wants CAP3
        
        // process options
    }
#endif
    
}
