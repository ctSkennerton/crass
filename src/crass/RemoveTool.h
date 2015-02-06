/*
 *  RemoveTool.h is part of the crisprtools project
 *  
 *  Created by Connor Skennerton on 22/12/11.
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

#ifndef crisprtools_RemoveTool_h
#define crisprtools_RemoveTool_h
#include <string>
#include <set>
#include <libcrispr/writer.h>

int removeMain(int argc, char ** argv);
void removeUsage(void);
int processRemoveOptions(int argc, char ** argv, std::set<std::string>& groups, std::string& outputFile, bool& remove );
void removeAssociatedData(xercesc::DOMElement * groupElement, crispr::xml::writer& xmlParser);
void parseMetadata(xercesc::DOMElement * parentNode, crispr::xml::writer& xmlParser);




#endif
