/*
 *  DeBruijn.cpp is part of the CRisprASSembler project
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

#include "DeBruijn.h"
#include "SeqUtils.h"
#include <libcrispr/Exception.h>

Kmer::Kmer(std::string str, bool t ) : Str(str), Terminal(t)
{}

void Kmer::print(std::ostream& out) {
	out << Str <<" : "<< Terminal<<std::endl;
}

DeBruijnGraph::DeBruijnGraph( int l) : Length(l){}

DeBruijnGraph::DeBruijnGraph( int l, std::map<StringToken, std::string>& seqs) : Length(l){
    generateGraph(seqs);
}


void DeBruijnGraph::generateGraph(std::map<StringToken, std::string>& seqs) {

    std::map<StringToken, std::string>::iterator iter;
    for (iter = seqs.begin(); iter != seqs.end(); iter++) {
        char** kmers = NULL;
        int * kmer_offsets = NULL;
		std::cout << iter->second <<std::endl;
        int num_mers = consumeSequence(iter->second, &kmers, &kmer_offsets);
        
        // add the kmers into the graph
        StringToken previous_token = 0;
        StringToken rc_previous_token = 0;
        std::string kmer_thread = "";
        std::string rc_kmer_thread = "";
        for(int i = 0; i < num_mers; ++i)
        {
            // make it a string!
            kmers[i][Length] = '\0';
            // check if the kmer is known
			StringToken token = IdConverter.getToken(kmers[i]);
            if (0 == token) {
                StringToken t = IdConverter.addString(kmers[i]);
                Kmer k = Kmer(kmers[i]);
                kmer_thread += to_string(t);
                
                if (i == 0 || i == num_mers - 1) {
                    k.setTerminal(true);
                }

                Nodes[t] = k;
                std::string rc = reverseComplement(kmers[i]);
                StringToken rct = IdConverter.addString(rc);
                Kmer rk = Kmer(rc);
                rc_kmer_thread += to_string(rct);
                if (i == 0 || i == num_mers - 1) {
                    rk.setTerminal(true);
                }
                Nodes[rct] = rk;

                if (previous_token != 0) {
                    Nodes[previous_token].addEdge(t);
                }
                if (rc_previous_token != 0) {
                    Nodes[rc_previous_token].addEdge(rct);
                }
                previous_token = t;
                rc_previous_token = rct;
            }
            else {
                if (previous_token != 0) {
                    Nodes[previous_token].addEdge(token);
                    previous_token = token;
                }
                if (rc_previous_token != 0) {
                    Nodes[rc_previous_token].addEdge(token);
                    rc_previous_token = token;
                }
            }
            //Nodes[previous_token].addDRType(iter->first, i);
        }
        Lookup[kmer_thread] = iter->second;
        Lookup[rc_kmer_thread] = iter->second;
        // clean up
        delete [] kmer_offsets;
        for(int i = 0; i < num_mers; i++)
        {
            delete [] kmers[i];
        }
        delete [] kmers;
        
    }
}


int DeBruijnGraph::consumeSequence(std::string& seq, char ***kmers, int **kmer_offsets){
    
    int str_len = static_cast<int>(seq.length());
    int off = str_len - Length;
    int num_mers = off + 1;
    
    // STOLED FROM SaSSY!!!!
    // First we cut kmers from the sequence then we use these to
    // determine overlaps, finally we make edges
    //
    // When we cut kmers from a read it is like this...
    //
    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    // ------------------------------
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    //
    // so we break the job into three parts
    //
    // XXXX|XXXXXXXXXXXXXXXXXXXXXX|XXXX
    // ----|----------------------|----
    // XXXX|XXXXXXXXXXXXXXXXXXXXXX|
    //  XXX|XXXXXXXXXXXXXXXXXXXXXX|X
    //   XX|XXXXXXXXXXXXXXXXXXXXXX|XX
    //    X|XXXXXXXXXXXXXXXXXXXXXX|XXX
    //     |XXXXXXXXXXXXXXXXXXXXXX|XXXX
    //
    // the first and last part may be a little slow but the middle part can fly through...
    //
    
    // make a 2d array for the kmers!
	//char ** kmers = NULL;
	//int * kmer_offsets = NULL;
	try {
		*kmers = new char*[num_mers];
	} catch(std::exception& e) {
		std::cerr<<"Attempting to alloc "<<num_mers<<std::endl;
		throw crispr::exception(__FILE__,
		                        __LINE__,
		                        __PRETTY_FUNCTION__,
		                        e.what());
	}
	try {
		for(int i = 0; i < num_mers; i++)
		{
			(*kmers)[i] = new char [Length+1];
		}
        // use these offsets when we cut kmers, they are a component of the algorithm
		(*kmer_offsets) = new int[num_mers];
		for(int i = 0; i < num_mers; i++)
		{
			(*kmer_offsets)[i] = i * -1; // Starts at [0, -1, -2, -3, -4, ...]
		}
	} catch(std::exception& e) {
		std::cerr<<"Attempting to alloc "<<Length+1<<std::endl;
		throw crispr::exception(__FILE__,
		                        __LINE__,
		                        __PRETTY_FUNCTION__,
		                        e.what());
	}
    int pos_counter = 0;
    
    // a slow-ish first part
    while(pos_counter < Length)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if(pos_counter >= j)
            {
                (*kmers)[j][(*kmer_offsets)[j]] = seq[pos_counter];
            }
            (*kmer_offsets)[j]++;
        }
        pos_counter++;
    }
    
    // this is the fast part of the loop
    while(pos_counter < off)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if((*kmer_offsets)[j] >= 0 && (*kmer_offsets)[j] < Length)
            {
                (*kmers)[j][(*kmer_offsets)[j]] = seq[pos_counter];
            }
            (*kmer_offsets)[j]++;
        }
        pos_counter++;
    }
    
    // an even slower ending
    while(pos_counter < str_len)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if((*kmer_offsets)[j] < Length)
            {
                (*kmers)[j][(*kmer_offsets)[j]] = seq[pos_counter];
            }
            (*kmer_offsets)[j]++;
        }
        pos_counter++;
    }

    return num_mers;
}

DeBruijnGraph::DataFound DeBruijnGraph::search(std::string seq){
    char** kmers = NULL;
    int * kmer_offsets = NULL;
    int num_mers = consumeSequence(seq, &kmers, &kmer_offsets);
    
    std::vector<StringToken> found_tokens;
    std::string token_thread = "";
    DataFound ret;
    ret.iFoundPosition = -1;
    //std::cout<< "--------" <<std::endl;
	for(int i = 0; i < num_mers; ++i)
    {
        // make it a string!
        kmers[i][Length] = '\0';
        StringToken token = IdConverter.getToken(kmers[i]);
        if(token != 0 && Nodes[token].getTerminal()){
            //found_tokens.push_back(token);
            token_thread += to_string(token);
			std::cout << "adding: "<< token<<std::endl;
//            if (Nodes[token].getLast()) {
//                // DR equal to the kmer length
//                std::cout << "DR in single kmer"<< std::endl;
//				ret.iFoundPosition = i;
//                ret.sDataFound = reconstructSequence(found_tokens);
//                return ret;            
//			}
            if (i != num_mers - 1) {
                int current_start = i;
                do {
                    ++i;
                    Kmer previous_kmer = Nodes[token];
                    token = IdConverter.getToken(kmers[i]);
                    if ((token != 0) & previous_kmer.hasEdge(token)) {
                        //found_tokens.push_back(token);
                        token_thread += to_string(token);
						std::cout << "adding: "<< token<<std::endl;
                        if (Nodes[token].getTerminal()) {
                            // reached the end of a DR
                            
                            // check to se if this thread is seen in the lookup
                            if (Lookup.find(token_thread) != Lookup.end()) {
                                ret.iFoundPosition = current_start;
                                ret.sDataFound = Lookup[token_thread];
                                return ret;
                            }                            
                        }
                    }
                    else {
                        // this kmer isn't in the graph
						std::cout<< "clearing kmer thread" <<std::endl;
                        //found_tokens.clear();
                        token_thread = "";
                        break;
                    }
                } while (i < num_mers);
            }
        }
    }
    return ret;
}

std::string DeBruijnGraph::reconstructSequence(std::vector<StringToken>& foundTokens){
    std::cout << "reconstructing sequence"<<std::endl;
	std::vector<StringToken>::iterator iter = foundTokens.begin();
    std::string reconstructed_sequence = Nodes[*iter].getSeq();
	std::cout << "creating first seed: "<< reconstructed_sequence<<std::endl;
    Nodes[*iter].print(std::cout);
    for (iter++; iter != foundTokens.end(); iter++) {
        std::string tmp = Nodes[*iter].getSeq();
		std::cout << "adding final base: " <<tmp[Length - 1] << " from kmer: "<< tmp<<std::endl;
        Nodes[*iter].print(std::cout);
        reconstructed_sequence += tmp[Length - 1];
    }
	std::cout << "final reconstructed sequence: "<< reconstructed_sequence<<std::endl;
    return reconstructed_sequence;
}

