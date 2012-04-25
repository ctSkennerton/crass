/*
 *  SearchChecker.h is part of the CRisprASSembler project
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

#ifndef crass_SearchChecker_h
#define crass_SearchChecker_h

#include <set>
#include <string>
#include <vector>
#include "ReadHolder.h"
#include "StringCheck.h"


class SearchData {
    ReadHolder * SD_ReadHolder;
    StringToken SD_Token;
    
    
public:
    SearchData() {
        SD_ReadHolder = NULL;
        SD_Token = 0;
    }
    
    inline StringToken token(){ return SD_Token;}
    inline void token(StringToken i) {SD_Token = i;}
    
    inline ReadHolder * read() {return SD_ReadHolder;}
    inline void read(ReadHolder * r) {SD_ReadHolder = r;}
};


class SearchChecker
{
    static SearchChecker *SC_instance;
    SearchChecker(){};
    std::string SC_FileName;             //a file containing a list of headers for interesting reads
    std::set<std::string> SC_Headers;    // a lookup for our interesting headers
    std::vector<SearchData> SC_Data;     // data about individual reads
    
public:
    // get/set
    inline bool hasHeader(std::string s) {return (SC_Headers.find(s) != SC_Headers.end());}
    inline void headerFile(std::string s ) {SC_FileName = s;}
    
    void processHeaderFile(void);
    static SearchChecker * instance();
};
static SearchChecker * debugger = SearchChecker::instance();

#endif
