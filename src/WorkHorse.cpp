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
#include "StringCheck.h"

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
        if (read_iter->second != NULL)
        {
            clearReadList(read_iter->second);
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
    std::vector<std::string>::iterator seq_iter = seqFiles.begin();
    while(seq_iter != seqFiles.end())
    {
        logInfo("Parsing file: " << *seq_iter, 1);
        
        // Need to make the string into something more old-skool so that
        // the search functions don't cry!
        char input_fastq[CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
        strncpy(input_fastq, seq_iter->c_str(), CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH);
        
        // direct repeat sequence and unique ID
        lookupTable patterns_lookup;
        
        // the sequence of whole spacers and their unique ID
        lookupTable reads_found;
        
        // Use a different search routine, depending on if we allow mismatches or not.
        READ_TYPE rt = decideWhichSearch(input_fastq);
        
        if(LONG_READ == rt)
        {
            mAveReadLength = crtSearchFastqFile(input_fastq, *mOpts, &mReads, &mStringCheck);
        }
        else
        {
            mAveReadLength = bmSearchFastqFile(input_fastq, *mOpts, patterns_lookup, reads_found, &mReads, &mStringCheck);
        }

        logInfo("Average read length: "<<mAveReadLength, 2);
        if (rt == LONG_READ) 
        {
            logInfo("long read algorithm selected", 2);
        }
        // only nessessary in instances where there are short reads
        if(SHORT_READ == rt)
        {
            findSingletons(input_fastq, *mOpts, patterns_lookup, reads_found, &mReads, &mStringCheck);
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

//**************************************
// Functions used to cluster DRs into groups and identify the "true" DR
//**************************************
int WorkHorse::mungeDRs(void)
{
    //-----
    // Cluster potential DRs and work out their true sequences
    // make the node managers while we're at it!
    //
    int next_free_GID = 1;
    std::map<std::string, int> k2GID_map;
    std::map<int, bool> groups;
    DR_Cluster_Map DR2GID_map;
    std::map<int, std::string> true_DRs;
    logInfo("Reducing list of potential DRs (1): Initial clustering", 1);
    logInfo("Reticulating splines...", 1);
    
    std::map<int, std::map<std::string, int> * > group_kmer_counts_map;

    // go through all of the read holder objects
    ReadMapIterator read_map_iter = mReads.begin();
    while (read_map_iter != mReads.end()) 
    {
        clusterDRReads(read_map_iter->first, &next_free_GID, &k2GID_map, &DR2GID_map, &groups, &group_kmer_counts_map);
        ++read_map_iter;
    }
    
    logInfo("Created: " << DR2GID_map.size() << " clusters", 2);
    logInfo("-------------", 4);
    DR_Cluster_MapIterator dcg_iter = DR2GID_map.begin();
    while(dcg_iter != DR2GID_map.end())
    {
        DR_ClusterIterator dc_iter = (dcg_iter->second)->begin();
        logInfo("Group: " << dcg_iter->first, 4);
        while(dc_iter != (dcg_iter->second)->end())
        {
            logInfo(mStringCheck.getString(*dc_iter), 4);
            dc_iter++;
        }
        dcg_iter++;
        logInfo("-------------", 4);
    }
    
    // print the reads to a file if requested
    if (mOpts->detect)
    {
        dumpReads(&DR2GID_map);
        return 1;
    }
    
    logInfo("Reducing list of potential DRs (2): Cluster refinement and true DR finding", 1);
    
    // go through all the counts for each group
    std::map<int, std::map<std::string, int> * >::iterator group_count_iter = group_kmer_counts_map.begin();
    while(group_count_iter != group_kmer_counts_map.end())
    {
        DR_Cluster * clustered_DRs = DR2GID_map[group_count_iter->first];
        if(clustered_DRs != NULL)
        {
            // it's real, so parse this group
            // get the five top kmers
            std::vector<std::string> n_top_kmers = getNMostAbundantKmers(CRASS_DEF_NUM_KMERS_4_MODE, group_count_iter->second);
            
            // a return value of false indicates that this function has deleted clustered_DRs
            parseGroupedDRs(group_count_iter->first, &n_top_kmers, clustered_DRs, &DR2GID_map, &next_free_GID, &true_DRs);
        }

        // delete the kmer count lists cause we're finsihed with them now
        if(NULL != group_count_iter->second)
        {
            delete group_count_iter->second;
            group_count_iter->second = NULL;
        }
        group_count_iter++;
    }

    // go through the DR2GID_map and make all reads in each group into nodes

    DR_Cluster_MapIterator drg_iter = DR2GID_map.begin();
    while(drg_iter != DR2GID_map.end())
    {
        if(NULL != drg_iter->second)
        {
            //std::cout << "Group: " << (drg_iter->first) << std::endl;
            
            mDRs[true_DRs[drg_iter->first]] = new NodeManager(true_DRs[drg_iter->first]);
            DR_ClusterIterator drc_iter = (drg_iter->second)->begin();
            while(drc_iter != (drg_iter->second)->end())
            {
                //std::cout << *drc_iter << std::endl;
                // go through each read
                
                drc_iter++;
            }
            //std::cout << "===================================" << std::endl;
        }
        drg_iter++;
    }
    return 0;
}

bool WorkHorse::parseGroupedDRs(int GID, std::vector<std::string> * nTopKmers, DR_Cluster * clustered_DRs, DR_Cluster_Map * DR2GID_map, int * nextFreeGID, std::map<int, std::string> * trueDRs)
{
    //-----
    // Cluster refinement and possible splitting for a Group ID
    //
    
    logInfo("Parsing group: " << GID, 4);
    
//++++++++++++++++++++++++++++++++++++++++++++++++
// Find a Master DR for this group of DRs

    // to store our DR which has all 5 kmers
    StringToken master_DR_token = -1;
    std::string master_DR_sequence = "**unset**";
    int master_read_count = 0;
    
    // these are needed fo rthe call to is kmer present but we don't actually need to values!
    bool disp_rc;
    int disp_pos;
    
    // go through the DRs in this cluster, we'd like to find one which has all kmers in it...
    // moreover, we need to get the one with all 5 and the most reads
    DR_ClusterIterator dr_iter = clustered_DRs->begin();
    while (dr_iter != clustered_DRs->end()) 
    {
        std::string tmp_DR = mStringCheck.getString(*dr_iter);
        bool got_all_mode_mers = true;
        std::vector<std::string>::iterator n_top_iter = nTopKmers->begin();
        while(n_top_iter != nTopKmers->end())
        {
            if(!isKmerPresent(&disp_rc, &disp_pos, &(*n_top_iter), &tmp_DR))
            {
                got_all_mode_mers = false;
                break;
            }
            n_top_iter++;
        }
        
        // did this guy have all 5?
        if(got_all_mode_mers)
        {
            int tmp_count = (int)(mReads[*dr_iter])->size();
            if(tmp_count > master_read_count)
            {
                master_read_count = tmp_count;
                master_DR_sequence = mStringCheck.getString(*dr_iter);
                master_DR_token = *dr_iter;
            }
        }
        
        // otherwise keep searching
        dr_iter++;
    }

    if(master_DR_token == -1)
    {
        // probably a dud. throw it out
        // free the memory and clean up
        logInfo("Could not identify a master DR", 4);
        delete clustered_DRs;
        clustered_DRs = NULL;
        (*DR2GID_map)[GID] = NULL;
        return false;
    }
    
    logInfo("Identified: " << master_DR_sequence << " as a master potential DR", 4);
    
    // now we have the 5 most abundant kmers and one DR which contains them all
    // time to rock and rrrroll!

//++++++++++++++++++++++++++++++++++++++++++++++++
// Initialise variables we'll need
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

//++++++++++++++++++++++++++++++++++++++++++++++++
// Set up the master DR's array and insert this guy into the main array
    
    // The offset of the start position of each potential DR 
    // when compared to the "true DR"
    // we use this structure when we detect overcollaping
    std::map<StringToken, int> DR_offset_map;
    
    // we will compare the orientations and positions of all other guys to the master so we need a permanent store
    bool kmer_rcs_DR_master[CRASS_DEF_NUM_KMERS_4_MODE];
    int kmer_positions_DR_master[CRASS_DEF_NUM_KMERS_4_MODE];
    
    // look for the start and end of the DR zone
    bool first_run = true;          // we only need to do this once
    int dr_zone_start = -1;
    int dr_zone_end = -1;
    
    // just the positions in the DR fangs...
    kmer_positions_ARRAY[0] = (int)(array_len/3);
    isKmerPresent(kmer_rcs_DR, kmer_positions_DR, &((*nTopKmers)[0]), &master_DR_sequence);
    
    for(int i = 1; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
    {
        isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), &((*nTopKmers)[i]), &master_DR_sequence);
        kmer_positions_ARRAY[i] = kmer_positions_DR[i] - kmer_positions_DR[0] + kmer_positions_ARRAY[0];
    }
    
    // store the first results away as the master results
    for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
    {
        kmer_rcs_DR_master[i] = kmer_rcs_DR[i];
        kmer_positions_DR_master[i] = kmer_positions_DR[i];
    }

    // note the position of the master DR in the array
    DR_offset_map[master_DR_token] = (kmer_positions_ARRAY[0] - kmer_positions_DR[0]);
    
    ReadListIterator read_iter = mReads[master_DR_token]->begin();
    while (read_iter != mReads[master_DR_token]->end()) 
    {
        // don't care about partials
        if(((*read_iter)->RH_StartStops[1] - (*read_iter)->RH_StartStops[0]) == (master_DR_sequence.length() - 1))
        {
            // the start of the read is 
            int this_read_start_pos = kmer_positions_ARRAY[0] - (*read_iter)->RH_StartStops[0] - kmer_positions_DR[0] ;
            if(first_run)
            {
                dr_zone_start =  this_read_start_pos + (*read_iter)->RH_StartStops[0];
                dr_zone_end =  this_read_start_pos + (*read_iter)->RH_StartStops[1];
                first_run = false;
            }
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

//++++++++++++++++++++++++++++++++++++++++++++++++
// now go thru all the other DRs in this group and add them into
// the consensus array
    dr_iter = clustered_DRs->begin();
    while (dr_iter != clustered_DRs->end()) 
    {
        // get the string for this mofo
        std::string tmp_DR = mStringCheck.getString(*dr_iter);
        
        // we've already done the master DR
        if(master_DR_token != *dr_iter)
        {
            // set this guy to -1 for now
            DR_offset_map[*dr_iter] = -1;
            
            // this is a DR we have yet to add to the coverage array
            // First we need to find the positions of the kmers in this DR
            for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
            {
                isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), &((*nTopKmers)[i]), &tmp_DR);
            }
            
            // we need to have at least half of the mode k-mers present to continue
            // Where they are not found, the position gets set to -1
            int num_neg1 = 0;
            for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
            {
                if(-1 == kmer_positions_DR[i])
                    num_neg1++; 
            }
            
            if(num_neg1 < CRASS_DEF_NUM_KMERS_4_MODE_HALF)
            {
                // all is good! ...so far
                // now we can determine if the DR is the reverse complement...
                int num_agree = 0;
                int num_disagree = 0;
                for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                {
                    if(kmer_rcs_DR_master[i] == kmer_rcs_DR[i])
                        num_agree++;
                    else
                        num_disagree++;
                }
                // if they are equal, do nothing
                if(num_agree != num_disagree)
                {
                    if(num_agree < num_disagree)
                    {
                        // we need to reverse all the reads and the DR for these reads
                        ReadListIterator read_iter = mReads[*dr_iter]->begin();
                        while (read_iter != mReads[*dr_iter]->end()) 
                        {
                            (*read_iter)->reverseComplementSeq();
                            read_iter++;
                        }
                        
                        // fix the places where the DR is stored
                        tmp_DR = reverseComplement(tmp_DR);
                        StringToken st = mStringCheck.addString(tmp_DR);
                        mReads[st] = mReads[*dr_iter];
                        mReads[*dr_iter] = NULL;
                        *dr_iter = st;
                    }

                    for(int i = 0; i < CRASS_DEF_NUM_KMERS_4_MODE; i++)
                    {
                        isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), &((*nTopKmers)[i]), &tmp_DR);
                    }
                    
                    // now it's time to try add this guy to the array...
                    // we need to find a suitable kmer...
                    int positioning_kmer_index = 0;
                    bool found_kmer = false;
                    // first find the first non -1 entry, there must be at least CRASS_DEF_NUM_KMERS_4_MODE_HALF
                    while(positioning_kmer_index < CRASS_DEF_NUM_KMERS_4_MODE)
                    {
                        if(-1 != kmer_positions_DR[positioning_kmer_index])
                            break;
                        positioning_kmer_index++;
                    }
                    // start here, now we look to the differences between this kmer and the next kmer
                    // and make sure that that difference is upheld in the master
                    while(positioning_kmer_index < (CRASS_DEF_NUM_KMERS_4_MODE - 1))
                    {
                        if((kmer_positions_DR[positioning_kmer_index] - kmer_positions_DR[positioning_kmer_index+1]) == (kmer_positions_DR_master[positioning_kmer_index] - kmer_positions_DR_master[positioning_kmer_index+1]))
                        {
                            found_kmer = true;
                            break;
                        }
                        positioning_kmer_index++;
                    }

                    if(found_kmer)
                    {
                        // note the position of this DR in the array
                        DR_offset_map[*dr_iter] = (kmer_positions_ARRAY[positioning_kmer_index] - kmer_positions_DR[positioning_kmer_index]);
                        
                        // We need to check that at least CRASS_DEF_PERCENT_IN_ZONE_CUT_OFF percent of bases agree within the "Zone"
                        int this_DR_start_index = 0;
                        int zone_start_index = dr_zone_start;
                        int comparison_length = (int)tmp_DR.length();
                        
                        // we only need to compare "within" the zone
                        if(DR_offset_map[*dr_iter] < dr_zone_start)
                        {
                            this_DR_start_index = dr_zone_start - DR_offset_map[*dr_iter];
                        }
                        else if(DR_offset_map[*dr_iter] > dr_zone_start)
                        {
                            zone_start_index = DR_offset_map[*dr_iter];
                        }
                        
                        // work out the comparison length
                        int eff_zone_length = dr_zone_end - zone_start_index;
                        int eff_DR_length = (int)tmp_DR.length() - this_DR_start_index;
                        if(eff_zone_length < eff_DR_length)
                            comparison_length = eff_zone_length;
                        else
                            comparison_length = eff_DR_length;

                        char cons_char = 'X';
                        int comp_end = zone_start_index + comparison_length;
                        double agress_with_zone = 0;
                        double comp_len = 0;
                        while(zone_start_index < comp_end)
                        {
                            // work out the consensus at this position
                            int max_count = 0;
                            for(int i = 0; i < 4; i++)
                            {
                                if(coverage_array[i][zone_start_index] > max_count)
                                {
                                    max_count = coverage_array[i][zone_start_index];
                                    cons_char = alphabet[i];
                                }
                            }
                            
                            // see if this DR agress with it
                            if(tmp_DR[this_DR_start_index] == cons_char)
                                agress_with_zone++;
                            comp_len++;
                            zone_start_index++;
                            this_DR_start_index++;
                        }

                        agress_with_zone /= comp_len;
                        if(agress_with_zone >= CRASS_DEF_PERCENT_IN_ZONE_CUT_OFF)
                        {
                            // we need to correct for the fact that we may not be using the 0th kmer
                            int positional_offset = kmer_positions_DR_master[0] - kmer_positions_DR_master[positioning_kmer_index] + kmer_positions_ARRAY[positioning_kmer_index];
                            ReadListIterator read_iter = mReads[*dr_iter]->begin();
                            while (read_iter != mReads[*dr_iter]->end()) 
                            {
                                // don't care about partials
                                if(((*read_iter)->RH_StartStops[1] - (*read_iter)->RH_StartStops[0]) == (tmp_DR.length() - 1))
                                {
                                    // we need to find the first kmer which matches the mode.
                                    int this_read_start_pos = positional_offset - (*read_iter)->RH_StartStops[0] -kmer_positions_DR[0] ;
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
                        }
                    }
                    else
                    {
                        // TODO
                        // should probably kill this guy here. He couldn't be added to the array
                    }
                }
            }
        }
        dr_iter++;
    }
//++++++++++++++++++++++++++++++++++++++++++++++++
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
        // we need at least CRASS_DEF_MIN_READ_DEPTH reads to call a DR
        if(total_count > CRASS_DEF_MIN_READ_DEPTH)
            conservation_array[j] = (float)(max_count)/total_count;
        else
            conservation_array[j] = 0;
    }
    

    for(int i = dr_zone_start; i <= dr_zone_end; i++)
    {
        std::cout << conservation_array[i] << ", ";
    }
    std::cout  << std::endl;
    for(int j = 0; j < 4; j++)
    {
        std::cout << alphabet[j] << ", ";
        for(int i = dr_zone_start; i <= dr_zone_end; i++)
        {
            std::cout << coverage_array[j][i] << ", ";
        }
        std::cout << std::endl;
    }
    
    // now use this information to find the true direct repeat
    // first work to the left
    while(dr_zone_start > 0)
    {
        if(conservation_array[dr_zone_start - 1] >= CRASS_DEF_CONS_CUT_OFF)
            dr_zone_start--;
        else
            break;
    }
    
    // next work to the right
    while(dr_zone_end < array_len - 1)
    {
        if(conservation_array[dr_zone_end + 1] >= CRASS_DEF_CONS_CUT_OFF)
            dr_zone_end++;
        else
            break;
    }
    
    // finally, make the true DR and check for consistency
    std::string true_DR = "";
    
    // use these variables to identify and store possible
    // collapsed clusters
    int collapsed_pos = -1;
    std::map<char, int> collapsed_options;            // holds the chars we need to split on
    std::map<int, bool> refined_DR_ends;              // so we can update DR ends based on consensus
    
    for(int i = dr_zone_start; i <= dr_zone_end; i++)
    {
        collapsed_pos++;
        if(conservation_array[i] >= CRASS_DEF_CONS_CUT_OFF)
        {
            refined_DR_ends[i] = true;
            true_DR += consensus_array[i];
        }
        else
        {
            // possible collapsed cluster
            refined_DR_ends[i] = false;
            logInfo("-------------", 5);
            logInfo("Possible collapsed cluster at position: " << collapsed_pos << " (" << (dr_zone_start + collapsed_pos) << ")", 5);
            
            float total_count = coverage_array[0][i] + coverage_array[1][i] + coverage_array[2][i] + coverage_array[3][i];
            logInfo("Base:  Count:  Cov:",5);
            for(int k = 0; k < 4; k++)
            {
                logInfo("  " << alphabet[k] << "     " << coverage_array[k][i] << "      " << ((float)coverage_array[k][i]/total_count), 5);
                if((float)coverage_array[k][i]/total_count >= CRASS_DEF_COLLAPSED_THRESHOLD)
                {
                    // there's enough bases here to warrant further investigation
                    collapsed_options[alphabet[k]] = (collapsed_options.size() + *nextFreeGID);
                    (*nextFreeGID)++;
                }
            }
            
            // make sure we've got more than 1 option
            if(2 > collapsed_options.size())
            {
                collapsed_options.clear();
                logInfo("", 5);
                logInfo("   ...ignoring", 5);
            }
            else
            {
                // make the collapsed pos array specific and exit this loop
                collapsed_pos += dr_zone_start;
                i = dr_zone_end + 1;
            }
        }
    }

    // check to make sure that the DR is not just some random long RE
    if(true_DR.length() > mOpts->highDRsize)
    {
        // probably a dud. throw it out
        // free the memory and clean up
        delete clustered_DRs;
        clustered_DRs = NULL;
        (*DR2GID_map)[GID] = NULL;
        
        logInfo("Killed: " << true_DR << " cause it was too long", 1);
        return false;
    }
    
//++++++++++++++++++++++++++++++++++++++++++++++++
// print out the consensus array 

    if(0 == collapsed_options.size())
    {
        // update the DR start and ends
        int diffs = dr_zone_end - dr_zone_start + 1 - (int)true_DR.length();
        while(0 != diffs)
        {
            // we need to update the start or end
            if(!refined_DR_ends[dr_zone_end])
            {
                dr_zone_end--;
                diffs--;
            }
            if(!refined_DR_ends[dr_zone_start])
            {
                dr_zone_start++;
                diffs--;
            }
        }
        
        if(isLogging(3))
        {
            stringstream ss;
            ss << std::endl << "%, ";
            for(int i = dr_zone_start; i <= dr_zone_end; i++)
            {
                ss << conservation_array[i] << ", ";
            }

            for(int j = 0; j < 4; j++)
            {
                ss << std::endl;
                ss << alphabet[j] << ", ";
                for(int i = dr_zone_start; i <= dr_zone_end; i++)
                {
                    ss << coverage_array[j][i] << ", ";
                }
            }
            logInfo(ss.str(), 3);
        }
    }

//++++++++++++++++++++++++++++++++++++++++++++++++
// clean up the mess we made

    delete[] consensus_array;
    delete[] conservation_array;
    for(int i = 0; i < 4; i++)
    {
        delete[] coverage_array[i];
    }
    
//++++++++++++++++++++++++++++++++++++++++++++++++
// possibly split the DR group

    if(collapsed_options.size() > 0)
    {
        // We need to build a bit of new infrastructure.
        // assume we have K different DR alleles and N putative DRs
        // we need to build K new DR clusters
        logInfo("", 5);
        logInfo("Attempting to split the collapsed DR", 5);
        std::map<char, int> coll_char_to_GID_map;
        std::map<char, int>::iterator co_iter = collapsed_options.begin();
        while(co_iter != collapsed_options.end())
        {
            int group = (*nextFreeGID)++;
            (*DR2GID_map)[group] = new DR_Cluster;
            coll_char_to_GID_map[co_iter->first] = group;
            co_iter++;
        }
    
        dr_iter = clustered_DRs->begin();
        while (dr_iter != clustered_DRs->end()) 
        {
            std::string tmp_DR = mStringCheck.getString(*dr_iter);
            if(-1 != DR_offset_map[*dr_iter])
            {
                // check if the deciding character is within range of this DR
                if(DR_offset_map[*dr_iter] <= collapsed_pos && collapsed_pos < (DR_offset_map[*dr_iter] + tmp_DR.length()))
                {
                    // this is easy, we can compare based on this char only
                    char decision_char = tmp_DR[collapsed_pos - DR_offset_map[*dr_iter]];
                    ((*DR2GID_map)[ coll_char_to_GID_map[ decision_char ] ])->push_back(*dr_iter);
                }
                else
                {
                    // this is tricky, we need to completely break the group and re-cluster
                    // from the ground up based on reads
                    // get the offset from the start of the DR to the deciding char
                    // if it is negative, the dec char lies before the DR
                    // otherwise it lies after
                    int dec_diff = collapsed_pos - DR_offset_map[*dr_iter];
                    
                    // we're not guaranteed to see all forms. So we need to be careful here...
                    // First we go through just to count the forms
                    std::map<char, ReadList *> forms_map;
                    
                    ReadListIterator read_iter = mReads[*dr_iter]->begin();
                    while (read_iter != mReads[*dr_iter]->end()) 
                    {
                        StartStopListIterator ss_iter = ((*read_iter)->RH_StartStops).begin();
                        while(ss_iter != ((*read_iter)->RH_StartStops).end())
                        {
                            int within_read_dec_pos = *ss_iter + dec_diff;
                            if(within_read_dec_pos > 0 && within_read_dec_pos < ((*read_iter)->RH_Seq).length())
                            {
                                char decision_char = ((*read_iter)->RH_Seq)[within_read_dec_pos];
                                forms_map[decision_char] = NULL;
                                break;
                            }
                            ss_iter+=2;
                        }
                        read_iter++;
                    }
                    
                    // the size of forms_map tells us how many different types we actually saw.
                    switch(forms_map.size())
                    {
                        case 1:
                        {
                            // we can just reuse the existing ReadList!
                            // find out which group this bozo is in
                            read_iter = mReads[*dr_iter]->begin();
                            bool break_out = false;
                            while (read_iter != mReads[*dr_iter]->end()) 
                            {
                                StartStopListIterator ss_iter = ((*read_iter)->RH_StartStops).begin();
                                while(ss_iter != ((*read_iter)->RH_StartStops).end())
                                {
                                    int within_read_dec_pos = *ss_iter + dec_diff;
                                    if(within_read_dec_pos > 0 && within_read_dec_pos < ((*read_iter)->RH_Seq).length())
                                    {
                                        char decision_char = ((*read_iter)->RH_Seq)[within_read_dec_pos];
                                        ((*DR2GID_map)[ coll_char_to_GID_map[ decision_char ] ])->push_back(*dr_iter);
                                        break_out = true;
                                        break;
                                    }
                                    ss_iter+=2;
                                }
                                read_iter++;             
                                if(break_out)     
                                    break;                  
                            }
                            break;
                        }
                        case 0:
                        {
                            // Somehting is wrong!
                            logWarn("No reads fit the form: " << tmp_DR, 1);
                            if(NULL != mReads[*dr_iter])
                            {
                                clearReadList(mReads[*dr_iter]);
                                delete mReads[*dr_iter];
                                mReads[*dr_iter] = NULL;
                            }
                            break;
                        }
                        default:
                        {
                            // we need to make a couple of new readlists and nuke the old one.
                            
                            // first make the new readlists
                            std::map<char, ReadList *>::iterator fm_iter = forms_map.begin();
                            while(fm_iter != forms_map.end())
                            {
                                // make the readlist
                                logInfo(tmp_DR << " : " << fm_iter->first << " : " << coll_char_to_GID_map[ fm_iter->first ], 1);
                                StringToken st = mStringCheck.addString(tmp_DR);
                                mReads[st] = new ReadList();
                                
                                // make sure we know which readlist is which
                                forms_map[fm_iter->first] = mReads[st];
                                
                                // put the new dr_token into the right cluster
                                ((*DR2GID_map)[ coll_char_to_GID_map[ fm_iter->first ] ])->push_back(st);
                                
                                // next!
                                fm_iter++;
                            }
                            
                            // put the correct reads on the correct readlist
                            read_iter = mReads[*dr_iter]->begin();
                            while (read_iter != mReads[*dr_iter]->end()) 
                            {
                                StartStopListIterator ss_iter = ((*read_iter)->RH_StartStops).begin();
                                while(ss_iter != ((*read_iter)->RH_StartStops).end())
                                {
                                    int within_read_dec_pos = *ss_iter + dec_diff;
                                    if(within_read_dec_pos > 0 && within_read_dec_pos < ((*read_iter)->RH_Seq).length())
                                    {
                                        char decision_char = ((*read_iter)->RH_Seq)[within_read_dec_pos];

                                        // push this readholder onto the correct list
                                        (forms_map[decision_char])->push_back(*read_iter);
                                        
                                        // make the original pointer point to NULL so we don't delete twice
                                        *read_iter = NULL;
                                        
                                        break;
                                    }
                                    ss_iter+=2;
                                }
                                read_iter++;                                    
                            }
                            
                            // nuke the old readlist
                            if(NULL != mReads[*dr_iter])
                            {
                                clearReadList(mReads[*dr_iter]);
                                delete mReads[*dr_iter];
                                mReads[*dr_iter] = NULL;
                            }                                
                            
                            break;
                        }
                    }
                }
            }
            dr_iter++;
        }
        
        // time to delete the old clustered DRs and the group from the DR2GID_map
        delete clustered_DRs;
        clustered_DRs = NULL;
        (*DR2GID_map)[GID] = NULL;
        
        logInfo("Calling the parser recursively", 4);
        
        // call this baby recursively with the new clusters
        std::map<char, int>::iterator cc_iter = coll_char_to_GID_map.begin();
        while(cc_iter != coll_char_to_GID_map.end())
        {
            parseGroupedDRs(cc_iter->second, nTopKmers, (*DR2GID_map)[cc_iter->second], DR2GID_map, nextFreeGID, trueDRs);
            cc_iter++;
        }
    }
    else
    {

        logInfo("Found DR: " << true_DR, 2);
        (*trueDRs)[GID] = true_DR;
    
//++++++++++++++++++++++++++++++++++++++++++++++++
// repair all the startstops for each read in this group
//
// This function is recursive, so we'll only get here when we have found exactly one DR

        DR_ClusterIterator drc_iter = ((*DR2GID_map)[GID])->begin();
        while(drc_iter != ((*DR2GID_map)[GID])->end())
        {
            // go through each read
            ReadListIterator read_iter = mReads[*drc_iter]->begin();
            while (read_iter != mReads[*drc_iter]->end()) 
            {
                std::cout << true_DR << std::endl << std::endl;
                std::cout << (*read_iter)->RH_Seq << std::endl;
                std::cout << (*read_iter)->splitApart() << std::endl;
                (*read_iter)->updateStartStops((DR_offset_map[*drc_iter] - dr_zone_start), &true_DR, mOpts);
                std::cout << (*read_iter)->splitApart() << std::endl;
                std::cout << ".............." << std::endl;
                read_iter++;
            }
            drc_iter++;
        }
    }
    
    return true;
}

bool WorkHorse::isKmerPresent(bool * didRevComp, int * startPosition, const std::string *  kmer, const std::string *  DR)
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
    std::string top_kmer;    
    std::map<std::string, bool> top_kmer_map;
    
    int iterations = (kmer_CountMap->size() < num2Get ) ? (int)kmer_CountMap->size() : num2Get;
    
    for (int i = 1; i <= iterations; i++) 
    {
        std::map<std::string, int>::iterator map_iter = kmer_CountMap->begin();
        int max_count = 0;
        
        while (map_iter != kmer_CountMap->end()) 
        {
            if((map_iter->second > max_count) && (top_kmer_map.find(map_iter->first) == top_kmer_map.end()))
            {
                max_count = map_iter->second;
                top_kmer = map_iter->first;
            }
            map_iter++;
        }
        top_kmer_map[top_kmer] = true;
    }
    std::vector<std::string> top_kmers;
    std::map<std::string, bool>::iterator tkm_iter = top_kmer_map.begin();
    while(tkm_iter != top_kmer_map.end())
    {
        top_kmers.push_back(tkm_iter->first);
        tkm_iter++;
    }
    return top_kmers;
}

bool WorkHorse::clusterDRReads(StringToken DRToken, int * nextFreeGID, std::map<std::string, int> * k2GIDMap, DR_Cluster_Map * DR2GIDMap, std::map<int, bool> * groups, std::map<int, std::map<std::string, int> * > * groupKmerCountsMap)
{
    //-----
    // hash a DR!
    //

    std::string DR = mStringCheck.getString(DRToken);
    int str_len = (int)DR.length();
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
        (*DR2GIDMap)[group] = new DR_Cluster;
        
        // we need a new kmer counter for this group
        (*groupKmerCountsMap)[group] = new std::map<std::string, int>;
    }
    
    // we need to record the group for this mofo!
    (*DR2GIDMap)[group]->push_back(DRToken);
    
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


void WorkHorse::dumpReads(DR_Cluster_Map * DR2GID_map)
{
    //-----
    // Print the reads from one cluster to a file...
    //

    DR_Cluster_MapIterator dr_clust_iter = DR2GID_map->begin();
    while (dr_clust_iter != DR2GID_map->end()) 
    {
        std::ofstream reads_file;
        reads_file.open((mOutFileDir + "Cluster_"+ to_string(dr_clust_iter->first) ).c_str());
        DR_ClusterIterator dr_iter = dr_clust_iter->second->begin();
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

