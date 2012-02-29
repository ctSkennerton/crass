/*
 * ExtractTool.h
 *
 * Copyright (C) 2011 2012 - Connor Skennerton
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
#ifndef EXTRACTTOOL_H
#define EXTRACTTOOL_H

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <bitset>
#include "Xml.h"

class ExtractTool 
{
	public:
    enum ELEMENT_TYPE{REPEAT,SPACER,CONSENSUS,FLANKER};
		// constructor
    ExtractTool();

		// destructor
		~ExtractTool();

		// option processing
		void generateGroupsFromString( std::string groupString);
		int processOptions(int argc, char ** argv);

		// process the input
		int processInputFile(const char * inputFile);
        void parseWantedGroups(crispr::XML& xmlObj, xercesc::DOMElement * rootElement);
        void extractDataFromGroup(crispr::XML& xmlDoc, xercesc::DOMElement * currentGroup);
    void processData(crispr::XML& xmlDoc, xercesc::DOMElement * currentType, ELEMENT_TYPE wantedType, std::string gid, std::ostream& outStream);
	private:
        std::set<std::string> ET_Group;             // holds a comma separated list of groups that need to be extracted
        std::ofstream ET_RepeatStream;
        std::ofstream ET_FlankerStream;
        std::ofstream ET_SpacerStream;
        std::ofstream ET_OneStream;
        std::ofstream ET_GroupStream;
        std::string ET_OutputPrefix;
        std::string ET_OutputNamePrefix;
        std::bitset<7> ET_BitMask;
/*
        each bit is for a different option:
bit pos:  6543210
          0000000
          |||||||
          ||||||------- a subset of the groups is wanted
          |||||-------- print different types of data (spacers, dr, flanker) into separate files
          ||||--------- print results from a different group into separate files
          |||---------- extract flanking sequences
          ||----------- extract spacer sequences
          |------------ extract direct repeat sequences
          ------------- print coverage information if extracting spacers
 */
};


/* 
	Global functions
*/

int extractMain(int argc, char ** argv);

void extractUsage(void);




#endif