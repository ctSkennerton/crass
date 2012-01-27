/*
 * MergeTool.h
 *
 * Copyright (C) 2011 - Connor Skennerton
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MERGETOOL_H
#define MERGETOOL_H
#include <set>
#include <string>
#include "Xml.h"
class MergeTool {
    std::set<std::string> MT_GroupIds;
    bool MT_Sanitise;
    int MT_NextGroupID;
    std::string MT_OutFile;
    crispr::XML MT_MasterDOM;
    
    
public:
    MergeTool(void){
        MT_OutFile = "crisprtools_merged.crispr";
        MT_NextGroupID = 1;
        MT_Sanitise = false;
    }
    
    ~MergeTool(){}
    inline bool getSanitise(void){return MT_Sanitise;};
    inline int getNextGroupID(void){return MT_NextGroupID;};
    inline void incrementGroupID(void){MT_NextGroupID++;};
    inline std::string getFileName(void){return MT_OutFile;};
    int processInputFile(const char * inputFile);
    int processOptions(int argc, char ** argv);
    
    inline std::set<std::string>::iterator find(std::string s) {return MT_GroupIds.find(s);};
    inline std::set<std::string>::iterator begin(){return MT_GroupIds.begin();};    
    inline std::set<std::string>::iterator end(){return MT_GroupIds.end();};
    inline void insert(std::string s){MT_GroupIds.insert(s);};
};


int mergeMain(int argc, char ** argv);

void mergeUsage(void);
#endif