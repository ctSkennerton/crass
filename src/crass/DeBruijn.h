/*
 *  DeBruijn.h is part of the CRisprASSembler project
 *
 *  Created by Connor Skennerton.
 *  Copyright 2011-2013 Connor Skennerton & Michael Imelfort. All rights reserved.
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

#ifndef __crass__DeBruijn__
#define __crass__DeBruijn__

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <set>
#include <libcrispr/StlExt.h>
#include "StringCheck.h"
#include "Types.h"


class Kmer {
    // The Id of this kmer
    std::string Str;
    // A vector of stringTokens that are edges to this kmer
    std::set<StringToken> Edges;
    // is this kmer the first in any DR type?
    bool Terminal;
public:
	Kmer(){};
    Kmer(std::string str, bool terminal = true);
    ~Kmer(){};
    inline void setTerminal(bool b) {Terminal = b;}
    inline bool getTerminal(void){return Terminal;}
    inline std::string getSeq(void){return Str;}

    
    inline void addEdge(StringToken edge) {Edges.insert(edge);}
    inline bool hasEdge(StringToken t) {return Edges.find(t) != Edges.end();}
    void print(std::ostream& out);
};

class DeBruijnGraph {

    // length of the Kmers in this graph
    int Length;
    // a map of all of the kmers in the graph
    std::map<StringToken, Kmer> Nodes;
    // a map of kmer id threads to kmer ids
    std::map<std::string, std::string> Lookup;
    // two-way hash joining Kmer sequences and kmer IDs
    StringCheck IdConverter;    
    
public:
    typedef struct _DataFound
	{
		int	            iFoundPosition;
		std::string     sDataFound;
	} DataFound;
    
    DeBruijnGraph( int l);
    DeBruijnGraph( int l, std::map<StringToken, std::string>& seqs);
    ~DeBruijnGraph(){};
    void generateGraph(std::map<StringToken, std::string>& seqs);
    int consumeSequence(std::string& seq, char ***kmers, int **kmer_offsets);
    DataFound search(std::string seq);
    std::string reconstructSequence(std::vector<StringToken>& foundTokens);
};
#endif /* defined(__crass__DeBruijn__) */
