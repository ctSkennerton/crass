/*
 *  Utils.cpp is part of the crisprtools project
 *  
 *  Created by Connor Skennerton on 3/12/11.
 *  Copyright 2011 Connor Skennerton. All rights reserved. 
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

#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <sys/stat.h>
#include <errno.h>
#include "Utils.h"
#include "StlExt.h"
#include "Exception.h"


void recursiveMkdir(std::string dir) 
{
    std::string tmp;
    size_t pos = 0;
    while ( std::string::npos != (pos = dir.find('/',pos+1)) ) {
        tmp = dir.substr(0,pos);
        mkdir(tmp.c_str(), (S_IRWXU | S_IRWXG));
    }
    mkdir(dir.c_str(), (S_IRWXU | S_IRWXG));
}

bool fileOrString(const char * str) {
    struct stat file_stats;
    if (-1 == stat(str, &file_stats)) {
        // some sort of error occured but 
        // it might still be a file as input
        // with just a weird error
        switch (errno) {
            case ENOENT:
                // no such file or directory
                // it's a string
                return false;
                break;
                
            default:
                // it was a file but there is some other error
                throw crispr::input_exception(strerror(errno));
                break;
        }
    } else {
        return true;
    }
}

void parseFileForGroups(std::set<std::string>& groups, const char * filePath) {
    // read through a file of group numbers
    
    std::string line;
    std::fstream in_file;
    
    in_file.open(filePath);
    if ( in_file.good()) {
        while (in_file >> line) {
            groups.insert(line);
        }
    } else {
        // error
        std::string s = "cannot read file ";
        throw crispr::input_exception( (s + filePath).c_str());
    }
}

void generateGroupsFromString(std::string str, std::set<std::string>& groups) {
    split(str, groups, ",");
}