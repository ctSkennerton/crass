// File: WorkHorse.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// This class is responsible for "running" the algorithm
// 
// --------------------------------------------------------------------
//  Copyright  2011 Michael Imelfort and Connor Skennerton
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//

#ifndef WorkHorse_h
    #define WorkHorse_h

// system includes
#include <iostream>
#include <string>
#include <vector>
#include <map>

// local includes
#include "crass_defines.h"
#include "libcrispr.h"
#include "NodeManager.h"
#include "ReadHolder.h"

// typedefs
typedef std::map<std::string, NodeManager *> DR_List;
typedef std::map<std::string, NodeManager *>::iterator DR_ListIterator;

typedef std::map<int, std::vector<std::string> *>::iterator DR_ClusterIterator;
typedef std::map<int, std::vector<std::string> *> DR_Cluster;

bool sortDirectRepeatByLength( const std::string &a, const std::string &b);


class WorkHorse {
    public:
        WorkHorse (const options * opts, std::string outDir) 
        { 
            mOpts = opts; 
            mOutFileDir = outDir;
            mAveReadLength = 0;
        }
        ~WorkHorse();
        
        // do all the work!
        int doWork(std::vector<std::string> seqFiles);

        //**************************************
        // file IO
        //**************************************
        //void write_spacerID(direct_repeat &dr_match, kseq_t *seq);
        //void write_direct_repeatID(direct_repeat &dr_match, kseq_t *seq);
        void writeLookupToFile(string &outFileName, lookupTable &outLookup);
        
    private:
        
        void clearReadList(ReadList * tmp_list);
        void clearReadMap(ReadMap * tmp_map);
        
        //**************************************
        // functions used to cluster DRs into groups and identify the "true" DR
        //**************************************
        int mungeDRs(void);                         // cluster potential DRs and make node managers
        bool clusterDRReads(std::string DR, int * nextFreeGID, std::map<std::string, int> * k2GIDMap, DR_Cluster * DR2GIDMap, std::map<int, bool> * groups, std::map<int, std::map<std::string, int> * > * groupKmerCountsMap);  // cut kmers and hash
        bool isKmerPresent(bool * didRevComp, int * startPosition, const std::string * kmer, const std::string * DR);
        std::vector<std::string> getNMostAbundantKmers(int num2Get, std::map<std::string, int> * kmer_CountMap);
        bool parseGroupedDRs(int GID, std::vector<std::string> * nTopKmers, std::vector<std::string> * clustered_DRs, DR_Cluster * DR2GID_map, int * nextFreeGID);
        
        //**************************************
        // file IO
        //**************************************
        void printFileLookups(std::string fileName, lookupTable &kmerLookup , lookupTable &patternsLookup, lookupTable &spacerLookup);
        void dumpReads(DR_Cluster * DR2GID_map);
        
    // members
        DR_List mDRs;                               // list of nodemanagers, cannonical DRs, one nodemanager per direct repeat
        ReadMap mReads;                             // reads containing possible double DRs
        const options * mOpts;                      // search options
        std::string mOutFileDir;                    // where to spew text to
        float mAveReadLength;                       // the average seen read length
};

#endif //WorkHorse_h