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
#include <algorithm>

// local includes
#include "WorkHorse.h"
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "NodeManager.h"
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "SmithWaterman.h"
#include "StlExt.h"

bool sortDirectRepeatByLength( const std::string &a, const std::string &b)
{
    return a.length() > b.length();
}

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
    clearReadMap( &mReads);
}

void WorkHorse::clearReadList(ReadList * tmp_list)
{
    //-----
    // clear all the reads from the readlist
    //
    ReadListIterator read_iter = tmp_list->begin();
    while(read_iter != tmp_list->end())
    {
        if(*read_iter != NULL)
        {
            delete *read_iter;
            *read_iter = NULL;
        }
        read_iter++;
    }
    tmp_list->clear();
}

void WorkHorse::clearReadMap(ReadMap * tmp_map)
{
    //-----
    // clear all the reads from the readlist
    //
    ReadMapIterator read_iter = tmp_map->begin();
    while(read_iter != tmp_map->end())
    {
        clearReadList(read_iter->second);
        if (read_iter->second != NULL)
        {
            delete read_iter->second;
            read_iter->second = NULL;
        }
        read_iter++;
    }
    tmp_map->clear();
}

// do all the work!
int WorkHorse::doWork(std::vector<std::string> seqFiles)
{
    //-----
    // Taken from the old main function
    //
    if(mOpts->max_mismatches == 0)
    {   
        logInfo("Finding CRISPRs using the boyer-moore search algorithm", 1); 
    }
    else
    {   
        logInfo("Finding CRISPRs using the bitap algorithm", 1); 
    }
    
    std::vector<std::string>::iterator seq_iter = seqFiles.begin();
    while(seq_iter != seqFiles.end())
    {
        logInfo("Parsing file: " << *seq_iter, 1);
        
        // Need to make the string into something more old-skool so that
        // the search functions don't cry!
        char input_fastq[CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
        strncpy(input_fastq, seq_iter->c_str(), CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH);
        
        // direct repeat sequence and unique ID
        lookupTable patternsLookup;
        
        // the sequence of whole spacers and their unique ID
        lookupTable readsFound;
        
        // Use a different search routine, depending on if we allow mismatches or not.
        if(mOpts->max_mismatches == 0)
        {   mAveReadLength = bmSearchFastqFile(input_fastq, *mOpts, patternsLookup, readsFound,  &mReads); }
        else
        {   mAveReadLength = bitapSearchFastqFile(input_fastq, *mOpts, patternsLookup, readsFound, &mReads); }

        logInfo("Average read length: "<<mAveReadLength, 2);
        
        // only nessessary in instances where there are short reads
        if (mAveReadLength < CRASS_DEF_READ_LENGTH_CUTOFF)
        {
            logInfo("Beginning multipattern matcher", 1);
            scanForMultiMatches(input_fastq, *mOpts, patternsLookup, readsFound, &mReads);
        }
        
        // There will be an abundance of forms for each direct repeat.
        // We needs to do somes clustering! Then trim and concatenate the direct repeats
        mungeDRs();
        
        // we don't use these tables any more so why print them
        //printFileLookups(*seq_iter, kmerLookup, patternsLookup, spacerLookup);
        
        logInfo("Finished file: " << *seq_iter, 1);
        seq_iter++;
    }
    
    logInfo("all done!", 1);
    return 0;
}

int WorkHorse::mungeDRs(void)
{
    //-----
    // Cluster potential DRs and work out their true sequences
    // make the node managers while we're at it!
    //
    logInfo("Reducing list of potential DRs", 1);
    
    int next_free_GID = 1;
    std::map<std::string, int> k2GID_map;
    std::map<int, bool> groups;
    DR_Cluster DR2GID_map;
    logInfo("Reducing list of potential DRs stage1", 1);
    
    std::map<int, std::map<std::string, int> * > group_kmer_counts_map;

    // go through all of the read holder objects
    ReadMapIterator read_map_iter = mReads.begin();
    while (read_map_iter != mReads.end()) 
    {
        clusterDRReads(read_map_iter->first, &next_free_GID, &k2GID_map, &DR2GID_map, &groups, &group_kmer_counts_map);
        ++read_map_iter;
    }
    
    logInfo("Reducing list of potential DRs stage2", 1);
    std::vector<std::string> mostOccuringKmers;
    
    // print the reads to a file if requested
    if (mOpts->detect)
    {
        dumpReads(&DR2GID_map);
        return 1;
    }
    
    // go through all the counts for each group
    std::map<int, std::map<std::string, int> * >::iterator group_count_iter = group_kmer_counts_map.begin();
    while(group_count_iter != group_kmer_counts_map.end())
    {
        std::vector<std::string> * clustered_DRs = DR2GID_map[group_count_iter->first];
        if(clustered_DRs != NULL)
        {
            // it's real!
            
            // get the five top kmers
            std::vector<std::string> n_top_kmers = getNMostAbundantKmers(CRASS_DEF_NUM_KMERS_4_MODE, group_count_iter->second);
            
            // go through the DRs in this cluster, we'd like to find one which has all 5 kmers in it...
            // moreover, we need to get the one with all 5 and the most reads
            std::cout << group_count_iter->first << std::endl;
            std::vector<std::string>::iterator dr_iter = clustered_DRs->begin();
            
            // to store our DR which has all 5 kmers
            std::string master_DR = "**unset**";
            int master_read_count = 0;
            
            // these are needed fo rthe call to is kmer present but we don't actually need to values!
            bool disp_rc;
            int disp_pos;
            while (dr_iter != clustered_DRs->end()) 
            {
                bool got_all_mode_mers = true;
                std::vector<std::string>::iterator n_top_iter = n_top_kmers.begin();
                while(n_top_iter != n_top_kmers.end())
                {
                    if(!isKmerPresent(&disp_rc, &disp_pos, &(*n_top_iter), &(*dr_iter)))
                    {
                        got_all_mode_mers = false;
                        break;
                    }
                    n_top_iter++;
                }
                
                // did this guy have all 5?
                if(got_all_mode_mers)
                {
                    int tmp_count = mReads[*dr_iter]->size();
                    if(tmp_count > master_read_count)
                    {
                        master_read_count = tmp_count;
                        master_DR = *dr_iter;
                    }
                }
                
                // otherwise keep searching
                dr_iter++;
            }
            std::cout << "-------------" << std::endl;
            
            if(master_DR == "**unset**")
            {
                // probably a dud. throw it out
                delete clustered_DRs;
                clustered_DRs = NULL;
            }
            else
            {
                std::cout << master_DR << std::endl;
                std::vector<std::string>::iterator n_top_iter = n_top_kmers.begin();
                while(n_top_iter != n_top_kmers.end())
                {
                    std::cout << *n_top_iter << " : " << (*(group_count_iter->second))[*n_top_iter] << std::endl;
                    n_top_iter++;
                }
            
                // now we have the 5 most abundant kmers and one DR which contains them all
                // time to rock and rrrroll!
                
                // chars we luv!
                char alphabet[4] = {'A', 'C', 'G', 'T'};
                
                // first we need a 4 * (3 * RL)
                int ** coverage_array = new int*[4];
                
                // fill it up!
                int array_len = 3 * (int)mAveReadLength;
                for(int i = 0; i < 4; i++)
                {
                    int * tmp_array = new int[array_len];
                    
                    // intialise to zeros!
                    for(int j = 0; j < array_len; j++)
                    {
                        tmp_array[j] = 0;
                    }
                    coverage_array[i] = tmp_array;
                }
                
                // we need a consensus array
                char * consensus_array = new char[array_len];
                for(int j = 0; j < array_len; j++)
                {
                    consensus_array[j] = 'X';
                }

                // we need a diversity array
                float * conservation_array = new float[array_len];
                for(int j = 0; j < array_len; j++)
                {
                    conservation_array[j] = 0;
                }
                
                // first we need to place the master DR about 1/3 of the way down
                // and then we need to know the positions of the kmers in the DR and in the 
                // coverage array
                int kmer_positions_DR[CRASS_DEF_NUM_KMERS_4_MODE];
                bool kmer_rcs_DR[CRASS_DEF_NUM_KMERS_4_MODE];
                int kmer_positions_ARRAY[CRASS_DEF_NUM_KMERS_4_MODE];
                for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                {
                    kmer_positions_DR[i] = -1;
                    kmer_rcs_DR[i] = false;
                    kmer_positions_ARRAY[i] = -1;
                }
                
                // just the positions in the DR fangs...
                kmer_positions_ARRAY[0] = (int)(array_len/3);
                isKmerPresent(kmer_rcs_DR, kmer_positions_DR, &(n_top_kmers[0]), &master_DR);
                
                for(int i = 1; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                {
                    isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), &(n_top_kmers[i]), &master_DR);
                    kmer_positions_ARRAY[i] = kmer_positions_DR[i] - kmer_positions_DR[0] + kmer_positions_ARRAY[0];
                }
                
                for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                {
                    std::cout << i << " : "  <<  kmer_positions_DR[i] << " : " <<  kmer_rcs_DR[i] << " : " << kmer_positions_ARRAY[i] << std::endl;
                }
                
                ReadListIterator read_iter = mReads[master_DR]->begin();
                while (read_iter != mReads[master_DR]->end()) 
                {
                    // don't care about partials
                    if(((*read_iter)->RH_StartStops[1] - (*read_iter)->RH_StartStops[0]) == (master_DR.length() - 1))
                    {
                        // the start of the read is 
                        int this_read_start_pos = kmer_positions_ARRAY[0] - (*read_iter)->RH_StartStops[0] -kmer_positions_DR[0] ;
                        for(int i = 0; i < ((*read_iter)->RH_Seq).length(); i++)
                        {
                            int index = -1;
                            switch((*read_iter)->RH_Seq[i])
                            {
                                case 'A':
                                    index = 0;
                                    break;
                                case 'C':
                                    index = 1;
                                    break;
                                case 'G':
                                    index = 2;
                                    break;
                                case 'T':
                                    index = 3;
                                    break;
                            }
                            if(index >= 0)
                            {
                                coverage_array[index][i+this_read_start_pos]++;
                            }
                        }
                    }
                    read_iter++;
                }
                
                // now we need to go thru all the other DRs in this group...
                dr_iter = clustered_DRs->begin();
                while (dr_iter != clustered_DRs->end()) 
                {
                    // we've already done the master DR
                    if(master_DR != *dr_iter)
                    {
                        // this is a DR we have yet to add to the coverage array
                        // First we need to find the positions of the kmers in this DR
                        std::cout << *dr_iter << std::endl;
                        for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                        {
                            isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), &(n_top_kmers[i]), &(*dr_iter));
                            std::cout << i << " : " << kmer_rcs_DR[i] << " : " <<  kmer_positions_DR[i] << std::endl;
                        }
                        ReadListIterator read_iter = mReads[*dr_iter]->begin();
                        while (read_iter != mReads[*dr_iter]->end()) 
                        {
                            // don't care about partials
                            if(((*read_iter)->RH_StartStops[1] - (*read_iter)->RH_StartStops[0]) == (master_DR.length() - 1))
                            {
/*                                int this_read_start_pos = kmer_positions_ARRAY[0] - (*read_iter)->RH_StartStops[0] -kmer_positions_DR[0] ;
                                for(int i = 0; i < ((*read_iter)->RH_Seq).length(); i++)
                                {
                                    int index = -1;
                                    switch((*read_iter)->RH_Seq[i])
                                    {
                                        case 'A':
                                            index = 0;
                                            break;
                                        case 'C':
                                            index = 1;
                                            break;
                                        case 'G':
                                            index = 2;
                                            break;
                                        case 'T':
                                            index = 3;
                                            break;
                                    }
                                    if(index >= 0)
                                    {
                                        coverage_array[index][i+this_read_start_pos]++;
                                    }
                                }*/
                            }
                            read_iter++;
                        }
                    }
                    dr_iter++;
                }
                
                
                // calculate consensus and diversity
                // warning, A heavy!
                for(int j = 0; j < array_len; j++)
                {
                    int max_count = 0;
                    float total_count = 0;
                    for(int i = 0; i < 4; i++)
                    {
                        total_count += (float)(coverage_array[i][j]);
                        if(coverage_array[i][j] > max_count)
                        {
                            max_count = coverage_array[i][j];
                            consensus_array[j] = alphabet[i];
                        }
                    }
                    if(0 != total_count)
                        conservation_array[j] = (float)(max_count)/total_count;
                }
                
                // print the arrays to std::out
                for(int i = 0; i < 4; i++)
                {
                    for(int j = 0; j < array_len; j++)
                    {
                        std::cout << coverage_array[i][j] << ",";
                    }
                    std::cout << std::endl;
                }
                
                for(int j = 0; j < array_len; j++)
                {
                    std::cout << consensus_array[j];// << ",";
                }
                std::cout << std::endl;
                
                for(int j = 0; j < array_len; j++)
                {
                    std::cout << conservation_array[j] << ",";
                }
                std::cout << std::endl;
                
                // clean up the mess we made
                delete[] consensus_array;
                delete[] conservation_array;
                for(int i = 0; i < 4; i++)
                {
                    delete[] coverage_array[i];
                }
            }
        }
        
        // delete the kmer count lists cause we're finsihed with them now
        if(NULL != group_count_iter->second)
            delete group_count_iter->second;
        group_count_iter->second = NULL;
        group_count_iter++;
    }
}



bool WorkHorse::isKmerPresent(bool * didRevComp, int * startPosition, const std::string * kmer, const std::string * DR)
{
    //-----
    // Work out if a Kmer is present in a string and store positions etc...
    //
    size_t pos = DR->find(*kmer);
    if(pos == string::npos)
    {
        // try the reverse complement
        // rev compt the kmer, it's shorter!
        std::string tmp_kmer = reverseComplement(*kmer);
        pos = DR->find(tmp_kmer);
        if(pos != string::npos)
        {
            // found the kmer!
            *didRevComp = true;
            *startPosition = (int)pos;           
            return true;
        }
    }
    else
    {
        // found the kmer!
        *didRevComp = false;
        *startPosition = (int)pos;
        return true;
    }
    *startPosition = -1;
    return false;
}

std::vector<std::string> WorkHorse::getNMostAbundantKmers(int num2Get, std::map<std::string, int> * kmer_CountMap)
{
    int previous_max = 100000000;
    std::string top_kmer;    
    std::vector<std::string> top_kmers;
    
    int iterations = (kmer_CountMap->size() < num2Get ) ? (int)kmer_CountMap->size() : num2Get;
    
    for (int i = 1; i <= iterations; i++) 
    {
        std::map<std::string, int>::iterator map_iter = kmer_CountMap->begin();
        int max_count = 0;
        
        while (map_iter != kmer_CountMap->end()) 
        {
            if (map_iter->second > max_count && map_iter->second < previous_max)
            {
                max_count = map_iter->second;
                top_kmer = map_iter->first;
            }
            map_iter++;
        }
        
        previous_max = max_count;
        top_kmers.push_back(top_kmer);
    }
    return top_kmers;
}

bool WorkHorse::clusterDRReads(std::string DR, int * nextFreeGID, std::map<std::string, int> * k2GIDMap, DR_Cluster * DR2GIDMap, std::map<int, bool> * groups, std::map<int, std::map<std::string, int> * > * groupKmerCountsMap)
{
    //-----
    // hash a DR!
    //
    int str_len = DR.length();
    int off = str_len - CRASS_DEF_KMER_SIZE;
    int num_mers = off + 1;
    
    //***************************************
    //***************************************
    //***************************************
    //***************************************
    // LOOK AT ME!
    // 
    // Here we declare the minimum criteria for membership when clustering
    // this is not cool!
    int min_clust_membership_count = mOpts->kmer_size;
    // 
    //***************************************
    //***************************************
    //***************************************
    //***************************************
    
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
    char ** kmers = new char*[num_mers];
    for(int i = 0; i < num_mers; i++)
    {
        kmers[i] = new char [CRASS_DEF_KMER_SIZE+1];
    }
    
    int * kmer_offsets = new int[num_mers];              // use these offsets when we cut kmers, they are a component of the algorithm
    for(int i = 0; i < num_mers; i++)
    {
        kmer_offsets[i] = i * -1; // Starts at [0, -1, -2, -3, -4, ...]
    }

    int pos_counter = 0;

    // a slow-ish first part
    while(pos_counter < CRASS_DEF_KMER_SIZE)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if(pos_counter >= j)
            {
                kmers[j][kmer_offsets[j]] = DR[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    
    // this is the fast part of the loop
    while(pos_counter < off)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if(kmer_offsets[j] >= 0 && kmer_offsets[j] < CRASS_DEF_KMER_SIZE)
            {
                kmers[j][kmer_offsets[j]] = DR[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    
    // an even slower ending
    while(pos_counter < str_len)
    {
        for(int j = 0; j < num_mers; j++)
        {
            if(kmer_offsets[j] < CRASS_DEF_KMER_SIZE)
            {
                kmers[j][kmer_offsets[j]] = DR[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    
    //
    // Now the fun stuff begins:
    //
    std::vector<std::string> homeless_kmers;
    std::map<int, int> group_count;
    std::map<std::string, int> local_kmer_CountMap;

    int group = 0;
    for(int i = 0; i < num_mers; ++i)
    {
        // make it a string!
        kmers[i][CRASS_DEF_KMER_SIZE+1]='\0';
        std::string tmp_str(kmers[i]);
        tmp_str = laurenize(tmp_str);
        
        // see if this guy has been counted!
        if(local_kmer_CountMap.find(tmp_str) == local_kmer_CountMap.end())
            local_kmer_CountMap[tmp_str] = 1;
        else
            local_kmer_CountMap[tmp_str]++;
        
        // see if we've seen this kmer before
        std::map<std::string, int>::iterator k2g_iter = k2GIDMap->find(tmp_str);
        if(k2g_iter == k2GIDMap->end())
        {
            // first time we seen this one
            homeless_kmers.push_back(tmp_str);
        }
        else
        {
            // only do this if our guy doesn't belong to a group yet
            if(0 == group)
            {
                // this kmer belongs to a group!
                std::map<int, int>::iterator this_group_iter = group_count.find(k2g_iter->second);
                if(this_group_iter == group_count.end())
                {
                    group_count[k2g_iter->second] = 1;
                }
                else
                {
                    group_count[k2g_iter->second]++;
                    if(min_clust_membership_count <= group_count[k2g_iter->second])
                    {
                        // we have found a group for this mofo!
                        group = k2g_iter->second;
                    }
                }
            }
        }
    }
    
    if(0 == group)
    {
        // we need to make a new group for all the homeless kmers
        group = (*nextFreeGID)++;
        
        // we need to make a new entry in the group map
        (*groups)[group] = true;
        (*DR2GIDMap)[group] = new std::vector<std::string>;
        
        // we need a new kmer counter for this group
        (*groupKmerCountsMap)[group] = new std::map<std::string, int>;
    }
    
    // we need to record the group for this mofo!
    (*DR2GIDMap)[group]->push_back(DR);
    
    // we need to assign all homeless kmers to the group!
    std::vector<std::string>::iterator homeless_iter = homeless_kmers.begin();
    while(homeless_iter != homeless_kmers.end())
    {
        (*k2GIDMap)[*homeless_iter] = group;
        homeless_iter++;
    }
    
    // we need to fix up the group counts
    std::map<std::string, int>::iterator local_count_iter = local_kmer_CountMap.begin();
    while(local_count_iter != local_kmer_CountMap.end())
    {
        (*(*groupKmerCountsMap)[group])[local_count_iter->first] += local_count_iter->second;
        local_count_iter++;
    }

    // clean up
    delete [] kmer_offsets;
    for(int i = 0; i < num_mers; i++)
    {
        delete [] kmers[i];
    }
    delete [] kmers;
    
    
}

//**************************************
// file IO
//**************************************

void WorkHorse::printFileLookups(std::string fileName, lookupTable &kmerLookup , lookupTable &patternsLookup, lookupTable &spacerLookup)
{
    //-----
    // Print all the information from a single round
    //
    logInfo("Printing lookup tables from file: " << fileName << "to " << mOpts->output_fastq, 1);
    
    // Make file names
    std::string kmer_lookup_name = mOpts->output_fastq + CRASS_DEF_DEF_KMER_LOOKUP_EXT;
    std::string patterns_lookup_name = mOpts->output_fastq + CRASS_DEF_DEF_PATTERN_LOOKUP_EXT;
    std::string spacer_lookup_name = mOpts->output_fastq + CRASS_DEF_DEF_SPACER_LOOKUP_EXT;
    
    // Write!
    writeLookupToFile(kmer_lookup_name, kmerLookup);  
    writeLookupToFile(patterns_lookup_name, patternsLookup);
    writeLookupToFile(spacer_lookup_name, spacerLookup);
}


void WorkHorse::dumpReads(DR_Cluster * DR2GID_map)
{
    //-----
    // Print the reads from one cluster to a file...
    //
    DR_ClusterIterator dr_clust_iter = DR2GID_map->begin();
    while (dr_clust_iter != DR2GID_map->end()) 
    {
        std::ofstream reads_file;
        reads_file.open((mOpts->output_fastq + "Cluster_"+ to_string(dr_clust_iter->first) ).c_str());
        std::vector<std::string>::iterator dr_iter = dr_clust_iter->second->begin();
        while (dr_iter != dr_clust_iter->second->end()) 
        {
            ReadListIterator read_iter = mReads[*dr_iter]->begin();
            while (read_iter != mReads[*dr_iter]->end()) 
            {
                reads_file<<">"<<(*read_iter)->RH_Header<<std::endl;
                reads_file<<(*read_iter)->RH_Seq<<std::endl;
                read_iter++;
            }
            dr_iter++;
        }
        reads_file.close();
        dr_clust_iter++;
    }
}


void WorkHorse::writeLookupToFile(string &outFileName, lookupTable &outLookup)
{
    std::ofstream outLookupFile;
    outLookupFile.open(outFileName.c_str());
    
    lookupTable::iterator ter = outLookup.begin();
    while (ter != outLookup.end()) 
    {
        outLookupFile<<ter->first<<"\t"<<ter->second<<endl;
        
        ter++;
    }
    outLookupFile.close();
}

