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

#include <map>
#include <string>
#include <vector>
#include <stack>
#include "ReadHolder.h"
#include "StringCheck.h"
//#include "libcrispr.h"
//#include "NodeManager.h"


typedef std::vector<StringToken> SearchDataNodes;

class SearchData {
    crass::ReadHolder * SD_holder;
    StringToken SD_drtoken;
    StringToken SD_nmHeaderToken;
    //NodeManager * SD_nodeManager;
    SearchDataNodes SD_nodes;
    std::string SD_truedr;
    int SD_groupNumber;
    std::vector<std::string> SD_spacerStringList;

public:
    SearchData() {
        SD_holder = NULL;
        //SD_nodeManager = NULL;
        SD_drtoken = 0;
        SD_nmHeaderToken = 0;
        SD_groupNumber = 0;
    }
    ~SearchData(){
        if (SD_holder != NULL) {
            delete SD_holder;
        }
    }
    inline StringToken token(){ return SD_drtoken;}
    inline void token(StringToken i) {SD_drtoken = i;}
    
    inline StringToken nmtoken() {return SD_nmHeaderToken;}
    inline void nmtoken(StringToken s) {SD_nmHeaderToken = s;}
    
    inline crass::ReadHolder * read() {return SD_holder;}
    inline void read(crass::ReadHolder * r) {SD_holder = r;}
    
    inline std::string truedr() {return SD_truedr;}
    inline void truedr(std::string s) {SD_truedr = s;}

    inline void addSpacer(std::string s) {SD_spacerStringList.push_back(s);}
    inline std::vector<std::string>::iterator beginSp() {return SD_spacerStringList.begin();}
    inline std::vector<std::string>::iterator endSp() {return SD_spacerStringList.end();}
    
    inline int gid() {return SD_groupNumber;}
    inline void gid(int i) {SD_groupNumber = i;}
    
    inline void addNode(StringToken t) {SD_nodes.push_back(t);}
    
    
    inline SearchDataNodes::iterator begin() {return SD_nodes.begin();}
    inline SearchDataNodes::iterator end() {return SD_nodes.end();}
};


typedef std::map<std::string, SearchData> SearchCheckerList;

class SearchChecker
{
    static SearchChecker *SC_instance;
    SearchChecker(){};
    std::string SC_FileName;             //a file containing a list of headers for interesting reads
    SearchCheckerList SC_Data;     // data about individual reads
    
public:
    // get/set
    inline bool hasHeader(std::string s) {return (find(s) != end());}
    inline void headerFile(std::string s ) {SC_FileName = s;}
    inline void add(std::string a, SearchData& s) {SC_Data[a] = s;}
    inline SearchCheckerList::iterator find(std::string s) {return SC_Data.find(s);}
    inline SearchCheckerList::iterator begin() {return SC_Data.begin();}
    inline SearchCheckerList::iterator end() {return SC_Data.end();}

    void processHeaderFile(void);
    static SearchChecker * instance();
};
static SearchChecker * debugger = SearchChecker::instance();

#endif
