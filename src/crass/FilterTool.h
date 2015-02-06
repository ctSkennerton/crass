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

#include <libcrispr/reader.h>
#include <bitset>
#include <set>
#include <vector>

class FilterTool {
    int FT_Spacers;
    int FT_Repeats;
    int FT_Flank;
    int FT_contigs;
    int FT_Coverage;
    std::string FT_OutputFile;
	int countElements(xercesc::DOMElement * parentNode);
   public: 
    FilterTool() {
        FT_Spacers = 0;
        FT_Repeats = 0;
        FT_Flank = 0;
        FT_contigs = 0;
        FT_Coverage = 0;
    }

int processOptions(int argc, char ** argv);
int processInputFile(const char * inputFile);
bool parseGroup(xercesc::DOMElement * parentNode, crispr::xml::reader& xmlParser);
bool parseData(xercesc::DOMElement * parentNode, crispr::xml::reader& xmlParser, std::set<std::string>& spacersToRemove);
inline int parseDrs(xercesc::DOMElement * parentNode){return countElements(parentNode);}
    int parseSpacers(xercesc::DOMElement * parentNode, crispr::xml::reader& xmlParser, std::set<std::string>& spacersToRemove);
inline int parseFlankers(xercesc::DOMElement * parentNode){return countElements(parentNode);}
    void parseAssembly(xercesc::DOMElement * parentNode, crispr::xml::reader& xmlParser, std::set<std::string>& spacersToRemove); 
    void parseContig(xercesc::DOMElement * parentNode, 
                                 crispr::xml::reader& xmlParser, 
                     std::string& contigId, std::set<std::string>& spacersToRemove);
    void parseCSpacer(xercesc::DOMElement * parentNode, 
                                  crispr::xml::reader& xmlParser, 
                      std::string& contigId, std::set<std::string>& spacersToRemove);
    void parseLinkSpacers(xercesc::DOMElement * parentNode, 
                                      crispr::xml::reader& xmlParser, 
                                      std::string& contigId, std::set<std::string>& spacersToRemove);
};

int filterMain(int argc, char ** argv);

void filterUsage(void);