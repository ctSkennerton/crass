/*
 * ExtractTool.h
 *
 * Copyright (C) 2011 - Connor Skennerton
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
#ifndef EXTRACTTOOL_H
#define EXTRACTTOOL_H

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include "Xml.h"

class ExtractTool 
{
	public:
    enum ELEMENT_TYPE{REPEAT,SPACER,CONSENSUS,FLANKER};
		// constructor
		ExtractTool();

		// destructor
		~ExtractTool();

		
		// get
		inline bool getDirectRepeat(void) { return ET_DirectRepeat;};
		inline bool getSpacer(void) { return ET_Spacer;};
		inline bool getFlanker(void) { return ET_Flanker;};
		inline std::set<std::string> getGroup(void) { return ET_Group;};
		inline bool getSplitGroup(void) { return ET_SplitGroup; };
		inline bool getSplitType(void) { return ET_SplitType;};
		inline bool getSubset(void) { return ET_Subset; };

		// set

		inline void setDirectRepeat(bool d) { ET_DirectRepeat = d;};
		inline void setSpacer(bool s) { ET_Spacer = s;};
		inline void setFlanker(bool f) { ET_Flanker = f;};
		inline void setGroup(std::set<std::string> g) {ET_Group = g;};
		inline void setSplitGroup(bool x) { ET_SplitGroup = x;};
		inline void setSplitType(bool y) {ET_SplitType = y;};

		// option processing
		void generateGroupsFromString( std::string groupString);
		int processOptions(int argc, char ** argv);

		// process the input
		int processInputFile(const char * inputFile);
        void parseWantedGroups(crispr::XML& xmlObj, xercesc::DOMElement * rootElement);
        void extractDataFromGroup(crispr::XML& xmlDoc, xercesc::DOMElement * currentGroup);
    void processData(crispr::XML& xmlDoc, xercesc::DOMElement * currentType, ELEMENT_TYPE wantedType, std::string gid, std::ostream& outStream);
	private:
		bool ET_DirectRepeat;						// extract the direct repeat sequences
		bool ET_Spacer;								// extract the spacer sequences
		bool ET_Flanker;							// extract the flanking sequences
		std::set<std::string> ET_Group;             // holds a comma separated list of groups that need to be extracted
		bool ET_SplitGroup; 						// print results for each group to a separate file
		bool ET_SplitType;							// print different types of results into different files
		bool ET_Subset;								// are we doing all groups
        std::ofstream ET_SpacerStream;              // the current output stream for spacers
        std::ofstream ET_RepeatStream;
        std::ofstream ET_FlankerStream;
        std::ofstream ET_OneStream;
        std::ofstream ET_GroupStream;
    std::string ET_OutputPrefix;
};


/* 
	Global functions
*/

int extractMain(int argc, char ** argv);

void extractUsage(void);




#endif