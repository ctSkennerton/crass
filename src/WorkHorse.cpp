// File: WorkHorse.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Implementation of WorkHorse functions
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
// system includes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <zlib.h>  
#include <fstream>

// local includes
#include "WorkHorse.h"
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "NodeManager.h"
#include "ReadHolder.h"

WorkHorse::~WorkHorse()
{
//    //-----
//    // destructor
//    //

    // clean up all the NodeManagers
    DR_ListIterator dr_iter = mDRs.begin();
    while(dr_iter != mDRs.end())
    {
        if(NULL != dr_iter->second)
        {
            delete dr_iter->second;
            dr_iter->second = NULL;
        }
        dr_iter++;
    }
    mDRs.clear();
    
    // clear the reads!
    clearReadLists();
}

void WorkHorse::clearReadLists(void)
{
    //-----
    // clear all the reads from the readlist
    //
    ReadListIterator read_iter = mReads.begin();
    while(read_iter != mReads.end())
    {
        if(*read_iter != NULL)
        {
            delete *read_iter;
            *read_iter = NULL;
        }
        read_iter++;
    }
    mReads.clear();
}

// do all the work!
int WorkHorse::doWork(std::vector<std::string> seqFiles)
{
    //-----
    // Taken from the old main function
    //
    if(mOpts->max_mismatches == 0)
    {   logInfo("Finding CRISPRs using the boyer-moore search algorithm", 1); }
    else
    {   logInfo("Finding CRISPRs using the bitap algorithm", 1); }
    
    std::vector<std::string>::iterator seq_iter = seqFiles.begin();
    while(seq_iter != seqFiles.end())
    {
        logInfo("Parsing file: " << *seq_iter, 1);
        
        // Need to make the string into something more old-skool so that
        // the search functions don't cry!
        char input_fastq[DEF_MCD_FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
        strncpy(input_fastq, seq_iter->c_str(), DEF_MCD_FASTQ_FILENAME_MAX_LENGTH);
        
        // vector of strings for out direct repeats 
        vector<string> patterns;
        
        float aveReadLength;
        
        // direct repeat sequence and unique ID
        lookupTable patternsLookup;
        
        // the sequence of the kmers on either side of the direct repeat and a unique ID
        lookupTable kmerLookup;
        
        // the sequence of whole spacers and their unique ID
        lookupTable spacerLookup;
        
        // Use a different search routine, depending on if we allow mismatches or not.
        if(mOpts->max_mismatches == 0)
        {   aveReadLength = bm_search_fastq_file(input_fastq, *mOpts, patternsLookup, spacerLookup, kmerLookup); }
        else
        {   aveReadLength = bitap_search_fastq_file(input_fastq, *mOpts, patternsLookup, spacerLookup, kmerLookup); }

        logInfo("Average read length: "<<aveReadLength, 2);
        
        // only nessessary in instances where there are short reads
        if (aveReadLength < DEF_MCD_READ_LENGTH_CUTOFF)
        {
            logInfo("Beginning multipattern matcher", 1);
            map_to_vector(patternsLookup, patterns);
            
            read_for_multimatch(input_fastq, *mOpts, patterns, kmerLookup);
        }
    
/*
********************************
add reads like this:

ReadHolder tmp_holder = new ReadHolder();
mReads->push_back(tmp_holder);

delete function is already present in the destructor!

********************************
*/
    
        // Checking step add here
        std::string dr_seq = "bleg";
        NodeManager * tmp_manager = new NodeManager(dr_seq);
        mDRs[dr_seq] = tmp_manager;

        printFileLookups(*seq_iter, kmerLookup, patternsLookup, spacerLookup);
        
        logInfo("Finished file: " << *seq_iter, 1);
        seq_iter++;
    }
    
    logInfo("all done!", 1);
    return 0;
}

//**************************************
// file IO
//**************************************

void WorkHorse::printFileLookups(std::string fileName, lookupTable &kmerLookup , lookupTable &patternsLookup, lookupTable &spacerLookup)
{
    //-----
    // Print all the information from a single round
    //
    logInfo("Printing lookup tables from file: " << fileName << "to " << mOutFileDir, 1);
    
    // Make file names
    string kmer_lookup_name = mOutFileDir + MCD_DEF_KMER_LOOKUP_EXT;
    string patterns_lookup_name = mOutFileDir + MCD_DEF_PATTERN_LOOKUP_EXT;
    string spacer_lookup_name = mOutFileDir + MCD_DEF_SPACER_LOOKUP_EXT;
    
    // Write!
    writeLookupToFile(kmer_lookup_name, kmerLookup);  
    writeLookupToFile(patterns_lookup_name, patternsLookup);
    writeLookupToFile(spacer_lookup_name, spacerLookup);
}

void WorkHorse::writeLookupToFile(string &outFileName, lookupTable &outLookup)
{
    ofstream outLookupFile;
    outLookupFile.open(outFileName.c_str());
    
    lookupTable::iterator ter = outLookup.begin();
    while (ter != outLookup.end()) 
    {
        outLookupFile<<ter->first<<"\t"<<ter->second<<endl;
        
        ter++;
    }
    outLookupFile.close();
}

