/*
 * SanitiseTool.h
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

#include <string>
#include <map>
#include <sstream>

#include "Xml.h"

typedef std::map<std::string, std::string> conversionMap;

class SanitiseTool {
    bool ST_Repeats;
    bool ST_Spacers;
    bool ST_contigs;
    bool ST_Flank;
    std::string ST_OutputFile;
    conversionMap ST_RepeatMap;
    conversionMap ST_SpacerMap;
    conversionMap ST_FlankMap;
    conversionMap ST_ContigMap;
    
    int ST_NextGroup;
    int ST_NextSpacer;
    int ST_NextRepeat;
    int ST_NextFlanker;
    int ST_NextContig;
    
    
public:
    // constructor
    SanitiseTool()
    {
        ST_NextGroup = 1;
        ST_NextSpacer = 1;
        ST_NextRepeat = 1;
        ST_NextFlanker = 1;
        ST_NextContig = 1;
        
        ST_Repeats = false;
        ST_Spacers = false;
        ST_contigs = false;
        ST_Flank = false;
    }
    
    
    
    // get and set
    inline int getNextGroup(void){return ST_NextGroup;}
    inline int getNextSpacer(void){return ST_NextSpacer;}
    inline int getNextRepeat(void){return ST_NextRepeat;}
    inline int getNextFlanker(void){return ST_NextFlanker;}
    inline int getNextContig(void){return ST_NextContig;}
    
    inline std::string getNextGroupS(void)
    {
        std::stringstream ss;
        ss<<'G'<<ST_NextGroup;
        return ss.str();
    }
    inline std::string getNextSpacerS(void)
    {
        std::stringstream ss;
        ss<<"SP"<<ST_NextSpacer;
        return ss.str();
    }
    inline std::string getNextRepeatS(void)
    {
        std::stringstream ss;
        ss<<"DR"<<ST_NextRepeat;
        return ss.str();
    }
    inline std::string getNextFlankerS(void)
    {
        std::stringstream ss;
        ss<<'F'<<ST_NextFlanker;
        return ss.str();
    }
    inline std::string getNextContigS(void)
    {
        std::stringstream ss;
        ss<<'C'<< ST_NextContig;
        return ss.str();
    }
    
    inline void setNextGroup(int i){ST_NextGroup = i;}
    inline void setNextSpacer(int i){ST_NextSpacer = i;}
    inline void setNextRepeat(int i){ST_NextRepeat = i;}
    inline void setNextFlanker(int i){ST_NextFlanker = i;}
    inline void setNextContig(int i){ST_NextContig = i;}
    
    
    inline void incrementGroup(void){ST_NextGroup++;}
    inline void incrementSpacer(void){ST_NextSpacer++;}
    inline void incrementRepeat(void){ST_NextRepeat++;}
    inline void incrementFlanker(void){ST_NextFlanker++;}
    inline void incrementContig(void){ST_NextContig++;}

    
    
    
    
    int processOptions(int argc, char ** argv);
    int processInputFile(const char * inputFile);
    void parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseDrs(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    
    void parseAssembly(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseContig(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseCSpacer(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseLinkSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseLinkFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
};


int sanitiseMain(int argc, char ** argv);

void sanitiseUsage(void);