/*
 * FilterTool.h
 *
 * Copyright (C) 2011, 2012 - Connor Skennerton
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

#include "Xml.h"
#include <bitset>

class FilterTool {
    int FT_Spacers;
    int FT_Repeats;
    int FT_Flank;
    int FT_contigs;
    std::string FT_OutputFile;
   public: 
    FilterTool() {
        FT_Spacers = 0;
        FT_Repeats = 0;
        FT_Flank = 0;
        FT_contigs = 0;
    }

int processOptions(int argc, char ** argv);
int processInputFile(const char * inputFile);
bool parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
bool parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
int parseDrs(xercesc::DOMElement * parentNode);
int parseSpacers(xercesc::DOMElement * parentNode);
int parseFlankers(xercesc::DOMElement * parentNode);
};

int filterMain(int argc, char ** argv);

void filterUsage(void);