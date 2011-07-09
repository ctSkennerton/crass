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
        
        // return value of the search functions
        float aveReadLength;
        
        // direct repeat sequence and unique ID
        lookupTable patternsLookup;
        
        
        // the sequence of whole spacers and their unique ID
        lookupTable readsFound;
        
        // Use a different search routine, depending on if we allow mismatches or not.
        if(mOpts->max_mismatches == 0)
        {   aveReadLength = bmSearchFastqFile(input_fastq, *mOpts, patternsLookup, readsFound,  &mReads); }
        else
        {   aveReadLength = bitapSearchFastqFile(input_fastq, *mOpts, patternsLookup, readsFound, &mReads); }

        logInfo("Average read length: "<<aveReadLength, 2);
        
        // only nessessary in instances where there are short reads
        if (aveReadLength < CRASS_DEF_READ_LENGTH_CUTOFF)
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
    logInfo("Reducing list of potential DRs1", 1);
    
    // go through all of the read holder objects
    ReadMapIterator read_map_iter = mReads.begin();
    while (read_map_iter != mReads.end()) 
    {
        clusterDRReads(read_map_iter->first, &next_free_GID, &k2GID_map, &DR2GID_map, &groups);
        ++read_map_iter;
    }
    logInfo("Reducing list of potential DRs2", 1);
    
    dumpReads(&DR2GID_map);
    
    // now that we know what our groups are it's time to find the one true direct repeat
    //oneDRToRuleThemAll(&DR2GID_map);
    logInfo("Reducing list of potential DRs3", 1);
    
    DR_ClusterIterator DR2GID_iter = DR2GID_map.begin();
    while(DR2GID_iter != DR2GID_map.end())
    {
        std::cout << DR2GID_iter->first << " : ";
        std::vector<std::string>::iterator vec_it;
        vec_it = DR2GID_iter->second->begin();
        while (vec_it != DR2GID_iter->second->end()) 
        {
            std::cout<<*vec_it<<" + ";
            vec_it++;
        }
        std::cout << std::endl;
        DR2GID_iter++;
    }
    logInfo("Reducing list of potential DRs4", 1);
    
    // create a crispr node object and add into a node manager
    
    logInfo("Done!", 1);
}


std::string WorkHorse::threadToSmithWaterman(std::vector<std::string> *array)
{
    //-----
    // take in a cluster of direct repeats and smithwaterman them to discover
    // the hidden powers within -- or just the correct direct repeat
    // if the group is small its totes okay to do all vs all
    // but we don't want to kill the machine so if it's a large group
    // we select some of the longer sequences in the group and perform the
    // comparison on just those sequences
    int group_size = array->size() - 1;
    if (CRASS_DEF_MAX_CLUSTER_SIZE_FOR_SW < group_size)
        group_size = CRASS_DEF_MAX_CLUSTER_SIZE_FOR_SW;

    // a hash of the alignments from smith waterman and their frequencies
    std::map<std::string, int> alignmentHash;
    
    for (int i = 0; i < group_size; i++) 
    {
        for (int j = i + 1; j <= group_size; j++) 
        {
           
            stringPair align_concensus = smithWaterman(array->at(i), array->at(j));
            
            // if the alignment length is less than CRASS_DEF_MIN_SW_ALIGNMENT_RATIO% of the string length
            if (align_concensus.first.length() < array->at(j).length() * CRASS_DEF_MIN_SW_ALIGNMENT_RATIO)
            {
            std::string rev_j = reverseComplement(array->at(j));
                stringPair rev_comp_align_concensus = smithWaterman(array->at(i), rev_j);
                
                if (rev_comp_align_concensus.first.length() > array->at(j).length() * CRASS_DEF_MIN_SW_ALIGNMENT_RATIO)
                {
                    if(1)//(rev_comp_align_concensus.first != rev_j) && (rev_comp_align_concensus.second != rev_j))
                    {
                        std::cout << rev_comp_align_concensus.first << " : " << rev_comp_align_concensus.second << " : " << array->at(i) << " : " <<  reverseComplement(array->at(j)) <<std::endl;
                        if(mReads.find(rev_comp_align_concensus.first) != mReads.end())
                        {
                            addOrIncrement(alignmentHash, rev_comp_align_concensus.first);
                        }
                        if(mReads.find(rev_comp_align_concensus.second) != mReads.end())
                        {
                            addOrIncrement(alignmentHash, rev_comp_align_concensus.second);
                        }
                    }
                }
            }
            else
            {
                if(1)//(align_concensus.first != array->at(j)) && (align_concensus.second != array->at(j)))
                {
                    std::cout << align_concensus.first << " : " << align_concensus.second << " : " << array->at(i) << " : " <<  array->at(j) <<std::endl;
                    if(mReads.find(align_concensus.first) != mReads.end())
                    {
                        addOrIncrement(alignmentHash, align_concensus.first);
                    }
                    if(mReads.find(align_concensus.second) != mReads.end())
                    {
                        addOrIncrement(alignmentHash, align_concensus.second);
                    }
                }
            }
        }
    }
    
    // The dr with the highest number of alignments *should* be the one we're looking for
    int max_val = 0;
    std::string the_true_DR = "**unset**";
    std::map<std::string, int>::iterator cluster_iter = alignmentHash.begin();
    while (cluster_iter != alignmentHash.end()) 
    {
        // use our kewl scoring code
        if(2 < cluster_iter->second)
        {
            int score = scorePotentialDR(cluster_iter->first, cluster_iter->second);
            if (score > max_val) 
            {
                max_val = score;
                the_true_DR = cluster_iter->first;
            }
        }
        cluster_iter++;
    }
    
    return the_true_DR;
}

int WorkHorse::scorePotentialDR(std::string DR, int multiplier)
{
    //-----
    // Calculate a "score" for each DR based on how
    // many times we saw it during the smith-waterman stage of the clustering
    // and also how many reads contained "doubled" versions of it
    //
    // get the Readlist for this fella
    ReadList * scoring_list = mReads.find(DR)->second;
    
    ReadListIterator score_iter = scoring_list->begin();
    ReadListIterator score_last = scoring_list->end();
    int dubs_count = 0;
    while(score_iter != score_last)
    {
        if(( ( (*score_iter)->RH_StartStops ).size() > 2 ) && !((((*score_iter)->RH_StartStops).front() == 0) || (((*score_iter)->RH_StartStops).back() == (*score_iter)->RH_Seq.length())))
        {
            dubs_count++;
        }
        score_iter++;
    }
    int score = DR.length() / multiplier * dubs_count;
    std::cout << DR << " : " << DR.length() << " / " << multiplier << " * " << dubs_count << " = " << score << std::endl;
    return score;
}


// use the true DR to fix the start and stop of the reads
void WorkHorse::clenseClusters(std::vector<std::string> * DR_array, std::string theTrueDR)
{
    // loop through all of the direct repeats in the cluster
    std::vector<std::string>::iterator DR_array_iter = DR_array->begin();
    
    while (DR_array_iter != DR_array->end()) 
    {
        // loop through all of the reads in that cluster
        ReadList * curr_list = mReads[*DR_array_iter];
                
        ReadListIterator curr_iter = curr_list->begin();

        while (curr_iter != curr_list->end()) 
        {
            // loop through the start stop list
            size_t array_length = (*curr_iter)->RH_StartStops.size() - 1;
            
            for (int dr_start_index_in_array = 0; dr_start_index_in_array < array_length; dr_start_index_in_array = dr_start_index_in_array + 2)
            {
                int aStartAlign, aEndAlign, aStartSearch, aEndSearch;
                
                // skip if front partial
//                if ( (*curr_iter)->RH_StartStops.at(dr_start_index_in_array) != 0 && (*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1) != (*curr_iter)->RH_Seq.length()) 
//                {
                    // make sure that we don't go off the beginning
                    int new_start = (*curr_iter)->RH_StartStops.at(dr_start_index_in_array) - CRASS_DEF_SW_SEARCH_EXT;
                    aStartSearch = (new_start >= 0) ? new_start : 0; 
                    
                    int new_end = (*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1) + CRASS_DEF_SW_SEARCH_EXT;
                    aEndSearch = (new_end <= (*curr_iter)->RH_Seq.length()) ? new_end : (*curr_iter)->RH_Seq.length();
                    
                    // search for the true direct repeat in the read within a window
                    smithWaterman((*curr_iter)->RH_Seq, theTrueDR, &aStartAlign, &aEndAlign, aStartSearch, aEndSearch);
                    
                logInfo("old start index: "<<(*curr_iter)->RH_StartStops.at(dr_start_index_in_array)<<" old end index: "<<(*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1), 6);
                
                    (*curr_iter)->RH_StartStops.at(dr_start_index_in_array) = aStartAlign;
                    (*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1) = aEndAlign;
                
                logInfo("new start index: "<<(*curr_iter)->RH_StartStops.at(dr_start_index_in_array)<<" new end index: "<<(*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1), 6);
                logInfo("the DR of the new index: "<<(*curr_iter)->RH_Seq.substr((*curr_iter)->RH_StartStops.at(dr_start_index_in_array), (*curr_iter)->RH_StartStops.at(dr_start_index_in_array + 1) - (*curr_iter)->RH_StartStops.at(dr_start_index_in_array)), 6);
                    // calculate the difference between the alignment start and end and the current numbers
//                }

                
            }
            curr_iter++;
        }
        DR_array_iter++;
    }
}

void WorkHorse::oneDRToRuleThemAll(DR_Cluster * DR2GID_map)
{
    //-----
    // Each DR_Cluster contains multiple variants on the true DR
    // But which one is real?
    //
    DR_ClusterIterator DR2GID_iter = DR2GID_map->begin();
    std::string the_true_DR = "unset";
    
    while (DR2GID_iter != DR2GID_map->end()) 
    {
        std::cout << DR2GID_iter->first << ":" << std::endl;
        
        // sort the vector from the longest to the shortest
        std::sort(DR2GID_iter->second->begin(), DR2GID_iter->second->end(), sortDirectRepeatByLength);
        
        //the_true_DR = threadToSmithWaterman(DR2GID_iter->second);
        
        
        //logInfo("Clustering has revealed the one true direct repeat: "<< the_true_DR, 5);
        //std::cout << "Clustering has revealed the one true direct repeat: " << the_true_DR << std::endl;

        //clenseClusters(DR2GID_iter->second, the_true_DR);

        ++DR2GID_iter;
    }
}

bool WorkHorse::clusterDRReads(std::string DR, int * nextFreeGID, std::map<std::string, int> * k2GIDMap, DR_Cluster * DR2GIDMap, std::map<int, bool> * groups)
{
    //-----
    // hash a DR!
    //
    int hash_size = 6;      // cutting 6 mers
    int str_len = DR.length();
    int off = str_len - hash_size;
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
        kmers[i] = new char [hash_size+1];
    }
    
    int * kmer_offsets = new int[num_mers];              // use these offsets when we cut kmers, they are a component of the algorithm
    for(int i = 0; i < num_mers; i++)
    {
        kmer_offsets[i] = i * -1; // Starts at [0, -1, -2, -3, -4, ...]
    }

    int pos_counter = 0;

    // a slow-ish first part
    while(pos_counter < hash_size)
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
            if(kmer_offsets[j] >= 0 && kmer_offsets[j] < hash_size)
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
            if(kmer_offsets[j] < hash_size)
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
    int group = 0;
    for(int i = 0; i < num_mers; ++i)
    {
        // make it a string!
        kmers[i][hash_size+1]='\0';
        std::string tmp_str(kmers[i]);
        tmp_str = laurenize(tmp_str);
        
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
    logInfo("Printing lookup tables from file: " << fileName << "to " << mOutFileDir, 1);
    
    // Make file names
    std::string kmer_lookup_name = mOutFileDir + CRASS_DEF_DEF_KMER_LOOKUP_EXT;
    std::string patterns_lookup_name = mOutFileDir + CRASS_DEF_DEF_PATTERN_LOOKUP_EXT;
    std::string spacer_lookup_name = mOutFileDir + CRASS_DEF_DEF_SPACER_LOOKUP_EXT;
    
    // Write!
    writeLookupToFile(kmer_lookup_name, kmerLookup);  
    writeLookupToFile(patterns_lookup_name, patternsLookup);
    writeLookupToFile(spacer_lookup_name, spacerLookup);
}


void WorkHorse::dumpReads(DR_Cluster * DR2GID_map)
{
    
    DR_ClusterIterator dr_clust_iter = DR2GID_map->begin();
    while (dr_clust_iter != DR2GID_map->end()) 
    {
        std::ofstream reads_file;
        reads_file.open((mOutFileDir + "Cluster_"+ to_string(dr_clust_iter->first) ).c_str());
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

