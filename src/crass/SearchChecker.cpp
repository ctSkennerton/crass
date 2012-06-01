/*
 *  SearchChecker.cpp is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011, 2012 Connor Skennerton & Michael Imelfort. All rights reserved. 
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
#include "SearchChecker.h"
#include <libcrispr/Exception.h>

SearchChecker * SearchChecker::SC_instance = NULL;


SearchChecker * SearchChecker::instance()
{
    if (!SC_instance)
        SC_instance = new SearchChecker;
        return SC_instance;
}

void SearchChecker::processHeaderFile() {
    std::fstream in;
    in.open(SC_FileName.c_str());
    if (in.good()) {
        std::string line;
        while (in >> line) {
            SearchData s;
            SC_Data[line] = s; 
        }
    } else {
        throw crispr::runtime_exception(__FILE__, 
                                        __LINE__,
                                        __PRETTY_FUNCTION__,
                                        "could not open header file");
    }
}