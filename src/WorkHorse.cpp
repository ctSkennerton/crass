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
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
// local includes
#include "WorkHorse.h"
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "NodeManager.h"
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "SmithWaterman.h"
#include "StlExt.h"
#include "StringCheck.h"
#include "config.h"
#ifdef PERFORM_CRASS_ASSEMBLY
#include "CrassXML.h"
#endif


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
    
    // check to see none of these are floating around
    DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(drg_iter->second != NULL)
		{
			delete drg_iter->second; 
		}
		drg_iter++;
	}
    
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

int WorkHorse::numOfReads(void)
{
    int count = 0;
    ReadMapIterator read_iter = mReads.begin();
    while(read_iter != mReads.end())
    {
        if (read_iter->second != NULL)
        {
            count += (int)(read_iter->second)->size();
        }
        read_iter++;
    }
    return count;
}

// do all the work!
int WorkHorse::doWork(std::vector<std::string> seqFiles)
{
    //-----
    // wrapper for the various processes needed to assemble crisprs
    //
	
	logInfo("Parsing reads in " << (seqFiles.size()) << " files", 1);
	if(parseSeqFiles(seqFiles))
	{
		logError("FATAL ERROR: parseSeqFiles failed");
        return 2;
	}
	
    // build the spacer end graph
	if(buildGraph())
	{
		logError("FATAL ERROR: buildGraph failed");
        return 3;
	}
	
#if DEBUG
	if (!mOpts->noDebugGraph) // this option will only exist if DEBUG is set anyway
    {
        // print debug graphs
        if(renderDebugGraphs())
        {
            logError("FATAL ERROR: renderDebugGraphs failed");
            return 4;
        }
    }

#endif
	
	// clean each spacer end graph
	if(cleanGraph())
	{
        logError("FATAL ERROR: cleanGraph failed");
        return 5;
	}
    
	// make spacer graphs
	if(makeSpacerGraphs())
	{
        logError("FATAL ERROR: makeSpacerGraphs failed");
        return 50;
	}
	
	// clean spacer graphs
	if(cleanSpacerGraphs())
	{
        logError("FATAL ERROR: cleanSpacerGraphs failed");
        return 51;
	}
	
	// make contigs
	if(splitIntoContigs())
	{
        logError("FATAL ERROR: splitIntoContigs failed");
        return 6;
	}
    
    //remove NodeManagers with low numbers of spacers
    if (removeLowSpacerNodeManagers())
    {
        logError("FATAL ERROR: removeLowSpacerNodeManagers failed");
        return 7;
    }
	
    if (generateFlankers()) {
        logError("FATAL ERROR: generateFlankers failed");
        return 70;
    }
    
    // dump the spacers to file
	if(dumpSpacers())
	{
        logError("FATAL ERROR: dumpSpacers failed");
        return 8;
	}
    
    // print the reads to a file if requested
	if(dumpReads(&mDR2GIDMap, false))
	{
        logError("FATAL ERROR: dumpReads failed");
        return 9;
	}
	
#if DEBUG
	if (!mOpts->noDebugGraph) 
    {
        // print clean graphs
        if(renderDebugGraphs("Clean_"))
        {
            logError("FATAL ERROR: renderDebugGraphs failed");
            return 10;
        }
    }

#endif
    
	// print spacer graphs
	if(renderSpacerGraphs())
	{
        logError("FATAL ERROR: renderSpacerGraphs failed");
        return 11;
	}
	
	printXML();
	
    logInfo("all done!", 1);
	return 0;
}

int WorkHorse::parseSeqFiles(std::vector<std::string> seqFiles)
{
	//-----
	// Load data from files and search for DRs
	//
    std::vector<std::string>::iterator seq_iter = seqFiles.begin();
    
    // direct repeat sequence and unique ID
    lookupTable patterns_lookup;
    
    // the sequence of whole spacers and their unique ID
    lookupTable reads_found;
    while(seq_iter != seqFiles.end())
    {
        logInfo("Parsing file: " << *seq_iter, 1);
        
        // Need to make the string into something more old-skool so that
        // the search functions don't cry!
        char input_fastq[CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
        strncpy(input_fastq, seq_iter->c_str(), CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH);
        
        // Use a different search routine, depending on if we allow mismatches or not.
        READ_TYPE rt = decideWhichSearch(input_fastq, &mAveReadLength, *mOpts);
        
        // Die if the average read length falls below 2DR + SP
        if (mAveReadLength < (2*mOpts->lowDRsize) + mOpts->lowSpacerSize) 
        {
            logInfo("The average read length "<<mAveReadLength<<" is below the minimum threshold of "<<(2*mOpts->lowDRsize) + mOpts->lowSpacerSize, 1);
            return 1;
        }
        
        if(rt == LONG_READ)
        {
            logInfo("Long read algorithm selected", 2);
            longReadSearch(input_fastq, *mOpts, &mReads, &mStringCheck, patterns_lookup, reads_found);
            logInfo("Number of reads found: "<<this->numOfReads(), 2);
            
        }
        else
        {
            logInfo("Short read algorithm selected", 2);
            shortReadSearch(input_fastq, *mOpts, patterns_lookup, reads_found, &mReads, &mStringCheck);
            logInfo("number of reads found so far: "<<this->numOfReads(), 2);
            
        }
        // Check to see if we found anything, should return if we haven't
        if (patterns_lookup.empty()) 
        {
            logInfo("No direct repeat sequences were identified for file: "<<input_fastq, 1);
        }
        logInfo("Finished file: " << *seq_iter, 1);
        
        seq_iter++;
    }
    
    if (patterns_lookup.size() > 0) 
    {
        seq_iter = seqFiles.begin();
        while (seq_iter != seqFiles.end()) {
            
            logInfo("Parsing file: " << *seq_iter, 1);
            
            // Need to make the string into something more old-skool so that
            // the search functions don't cry!
            char input_fastq[CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH] = { '\0' };
            strncpy(input_fastq, seq_iter->c_str(), CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH);
            
            logInfo("Begining Second iteration through files to recruit singletons", 2);
            
            findSingletons(input_fastq, *mOpts, patterns_lookup, reads_found, &mReads, &mStringCheck);
            
            seq_iter++;
        }
    }
    logInfo("Searching complete. " << mReads.size()<<" direct repeat variants have been found", 1);
    logInfo("number of reads found so far: "<<this->numOfReads(), 2);
    
    // There will be an abundance of forms for each direct repeat.
    // We needs to do somes clustering! Then trim and concatenate the direct repeats
    if (mungeDRs())
    {
        logError("Wierd stuff happend when trying to get the 'true' direct repeat");            
    }
    
    return 0;
}

int WorkHorse::buildGraph(void)
{
	//-----
	// Load the spacers into a graph
	//
    // go through the DR2GID_map and make all reads in each group into nodes
    DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
    while(drg_iter != mDR2GIDMap.end())
    {
        if(NULL != drg_iter->second)
        {            
#ifdef DEBUG
            logInfo("Creating NodeManager "<<drg_iter->first, 6);
#endif
            mDRs[mTrueDRs[drg_iter->first]] = new NodeManager(mTrueDRs[drg_iter->first], mOpts);
            DR_ClusterIterator drc_iter = (drg_iter->second)->begin();
            while(drc_iter != (drg_iter->second)->end())
            {
                // go through each read
                ReadListIterator read_iter = mReads[*drc_iter]->begin();
                while (read_iter != mReads[*drc_iter]->end()) 
                {
                    mDRs[mTrueDRs[drg_iter->first]]->addReadHolder(*read_iter);
                    read_iter++;
                }
                drc_iter++;
            }
        }
        drg_iter++;
    }
    return 0;
}

int WorkHorse::cleanGraph(void)
{
	//-----
	// Wrapper for graph cleaning
	//
	logInfo("Cleaning graphs", 1);
	DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(NULL != drg_iter->second)
		{            
#ifdef DEBUG
            if (NULL == mDRs[mTrueDRs[drg_iter->first]])
            {
                logWarn("Before Clean Graph: NodeManager "<<drg_iter->first<<" is NULL",6);
            }
            else
            {
#endif
                if((mDRs[mTrueDRs[drg_iter->first]])->cleanGraph())
                {
                    return 1;
                }
#ifdef DEBUG
            }
            if (NULL == mDRs[mTrueDRs[drg_iter->first]])
            {
                logWarn("After Clean Graph: NodeManager "<<drg_iter->first<<" is NULL",6);
            }
#endif
		}
		drg_iter++;
	}
	return 0;
}

int WorkHorse::removeLowSpacerNodeManagers(void)
{
    logInfo("Removing CRISPRs with low numbers of spacers", 1);
	DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(NULL != drg_iter->second)
		{            
#ifdef DEBUG
            if (NULL == mDRs[mTrueDRs[drg_iter->first]])
            {
                logWarn("Before Low Spacer Removal: NodeManager "<<drg_iter->first<<" is NULL",6);
            }
            else
            {
#endif
                if((mDRs[mTrueDRs[drg_iter->first]])->getSpacerCount(false) < mOpts->covCutoff) 
                {
                    logInfo("Deleting NodeManager "<<drg_iter->first<<" as it contained less than "<<mOpts->covCutoff<<" attached spacers",4);
                    delete mDRs[mTrueDRs[drg_iter->first]];
                    mDRs[mTrueDRs[drg_iter->first]] = NULL;
                }
#ifdef DEBUG
            }
#endif
            
		}
		drg_iter++;
	}
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
    logInfo("Reducing list of potential DRs (1): Initial clustering", 1);
    logInfo("Reticulating splines...", 1);
    
    std::map<int, std::map<std::string, int> * > group_kmer_counts_map;
    
    // go through all of the read holder objects
    ReadMapIterator read_map_iter = mReads.begin();
    while (read_map_iter != mReads.end()) 
    {
        clusterDRReads(read_map_iter->first, &next_free_GID, &k2GID_map, &group_kmer_counts_map);
        ++read_map_iter;
    }
    
    logInfo("Reducing list of potential DRs (2): Purging singleton and low abundance clusters", 1);
    int purge_counter = 0;
    DR_Cluster_MapIterator dcg_iter = mDR2GIDMap.begin();
    while(dcg_iter != mDR2GIDMap.end())
    {
        // delete groups that only have a single direct repeat variant as they should be bad!
        if((dcg_iter->second)->size() == 1)
        {
            logInfo("Purging Group " << dcg_iter->first<<" as it contains a single DR variant", 4);
            // make this group point to null so that we won't go through it again
            if(NULL != dcg_iter->second)
            {
            	delete dcg_iter->second;
            	dcg_iter->second = NULL;
            }
            if(NULL != mDR2GIDMap[dcg_iter->first])
            {
            	delete mDR2GIDMap[dcg_iter->first];
            	mDR2GIDMap[dcg_iter->first] = NULL;
            }
            purge_counter++;
        }
        else if(isLogging(4))
        {
            DR_ClusterIterator dc_iter = (dcg_iter->second)->begin();
            if (dcg_iter->second != NULL) 
            {
                logInfo("-------------", 4);
                logInfo("Group: " << dcg_iter->first, 4);
                while(dc_iter != (dcg_iter->second)->end())
                {
                    logInfo(mStringCheck.getString(*dc_iter), 4);
                    dc_iter++;
                }
                logInfo("-------------", 4);
            }
            
        }
        dcg_iter++;
    }
    
    
    if (purge_counter == 1) 
    {
        logInfo(purge_counter<<" cluster was purged",2);
    } 
    else 
    {
        logInfo(purge_counter<<" clusters were purged",2);
        
    }
    
    logInfo("Reducing list of potential DRs (3): Cluster refinement and true DR finding", 1);
    
    // go through all the counts for each group
    std::map<int, std::map<std::string, int> * >::iterator group_count_iter = group_kmer_counts_map.begin();
    while(group_count_iter != group_kmer_counts_map.end())
    {
        DR_Cluster * clustered_DRs = mDR2GIDMap[group_count_iter->first];
        if(clustered_DRs != NULL)
        {
            // it's real, so parse this group
            // get the five top kmers
            std::vector<std::string> n_top_kmers;
            if (getNMostAbundantKmers(n_top_kmers, CRASS_DEF_NUM_KMERS_4_MODE, group_count_iter->second)) 
            {
                // a return value of false indicates that this function has deleted clustered_DRs
                parseGroupedDRs(group_count_iter->first, &n_top_kmers, clustered_DRs, &next_free_GID);
            }
            else
            {
                // if there is less kill the group
                if(NULL != group_count_iter->second)
                {
                	delete group_count_iter->second;
                	group_count_iter->second = NULL;
                }
                if(NULL != mDR2GIDMap[group_count_iter->first])
                {
                	delete mDR2GIDMap[group_count_iter->first];
                	mDR2GIDMap[group_count_iter->first] = NULL;
                }
            }
        }
        
        // delete the kmer count lists cause we're finsihed with them now
        if(NULL != group_count_iter->second)
        {
            delete group_count_iter->second;
            group_count_iter->second = NULL;
        }
        group_count_iter++;
    }
    
    return 0;
}

bool WorkHorse::parseGroupedDRs(int GID, std::vector<std::string> * nTopKmers, DR_Cluster * clustered_DRs, int * nextFreeGID)
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
        
        // did this guy have all n?
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
        if(NULL != clustered_DRs)
        {
        	delete clustered_DRs;
        	clustered_DRs = NULL;
            mDR2GIDMap.erase(GID);
        }
//        if(NULL != mDR2GIDMap[GID])
//        {
//        	//delete mDR2GIDMap[GID];
//        	//mDR2GIDMap[GID] = NULL;
//        }
        return false;
    }
    
    logInfo("Identified: " << master_DR_sequence << " as a master potential DR", 4);
    
    // now we have the n most abundant kmers and one DR which contains them all
    // time to rock and rrrroll!
    
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialise variables we'll need
    // chars we luv!
    char alphabet[4] = {'A', 'C', 'G', 'T'};
    
    // first we need a 4 * array_len
    int ** coverage_array = new int*[4];
    
    int array_len = ((int)3*mAveReadLength > 1200) ? (int)3*mAveReadLength : 1200;
    
    // fill it up!
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
        int dr_start_index = 0;
        int dr_end_index = 1;
        
        // If you -1 from the master_DR_sequence.length() this loop will through an error as it will never be true
        while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) != ((int)(master_DR_sequence.length()) - 1))
        {
            dr_start_index += 2;
            dr_end_index += 2;
        }
        // this should be OK
        // If you -1 from the master_DR_sequence.length() this loop will throw an error as it will never be true
        if(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == ((int)(master_DR_sequence.length()) - 1))
        {
            // the start of the read is 
        	//MI logInfo(dr_start_index << " : "  << kmer_positions_ARRAY[0] << " : " << (*read_iter)->getSeq() << " : " << (*read_iter)->getSeqCharAt((*read_iter)->startStopsAt(dr_start_index))<< " : " << (*read_iter)->startStopsAt(dr_start_index) << " : " << kmer_positions_DR[0], 1);
        	//MI logInfo((*read_iter)->splitApartSimple(),1);
        	//MI logInfo((*read_iter)->DRLowLexi(),1);
            int this_read_start_pos = kmer_positions_ARRAY[0] - (*read_iter)->startStopsAt(dr_start_index) - kmer_positions_DR[0] ;
            if(first_run)
            {
                dr_zone_start =  this_read_start_pos + (*read_iter)->startStopsAt(dr_start_index);
                dr_zone_end =  this_read_start_pos + (*read_iter)->startStopsAt(dr_end_index);
                first_run = false;
            }
            for(int i = 0; i < (int)(*read_iter)->getSeqLength(); i++)
            {
                int index = -1;
                switch((*read_iter)->getSeqCharAt(i))
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
                	if((i+this_read_start_pos) >= array_len)
                		logError("The consensus/coverage arrays are too short. Consider changing the array_len variable to something larger");
                    coverage_array[index][i+this_read_start_pos]++;
                }
            }
        }
        else
        {
            logError("Everything is wrong (A)");
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
                // there's not too much we can do
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
                                int dr_start_index = 0;
                                int dr_end_index = 1;
                                while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) != ((int)(tmp_DR.length()) - 1))
                                {
                                    dr_start_index += 2;
                                    dr_end_index += 2;
                                }        
                                do
                                {
									if(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == (((int)(tmp_DR.length())) - 1))
									{
										// we need to find the first kmer which matches the mode.
										int this_read_start_pos = positional_offset - (*read_iter)->startStopsAt(dr_start_index) - kmer_positions_DR[0] ;
										for(int i = 0; i < (int)(*read_iter)->getSeqLength(); i++)
										{
											int index = -1;
											switch((*read_iter)->getSeqCharAt(i))
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
							                	if((i+this_read_start_pos) >= 0)
							                	{
							                		coverage_array[index][i+this_read_start_pos]++;
							                	}
											}
										}
									}
                                    // go onto the next DR
                                    dr_start_index += 2;
                                    dr_end_index += 2;
                                    
                                    // check that this makes sense
                                    if(dr_start_index >= (int)((*read_iter)->numRepeats()*2))
                                    	break;
                                    
                                } while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == (((int)(tmp_DR.length())) - 1));
                                //}
                                //else
                                //  logError("Everything is wrong (B)");
                                read_iter++;
                            }
                        }
                    }
                    else
                    {
                        std::cout<<"TODO: kill this guy here"<<std::endl;
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
    
    // trim these back a bit (if we trim too much we'll get it back right now anywho)
    dr_zone_start += CRASS_DEF_DR_ZONE_TRIM_AMOUNT;
    dr_zone_end -= CRASS_DEF_DR_ZONE_TRIM_AMOUNT;
    
    // now use this information to find the true direct repeat
    // first work to the left
    while(dr_zone_start > 0)
    {
        if(conservation_array[dr_zone_start - 1] >= CRASS_DEF_ZONE_EXT_CONS_CUT_OFF)
        {
            dr_zone_start--;
        }
        else
        {
            break;
        }
    }
    
    // next work to the right
    while(dr_zone_end < array_len - 1)
    {
        if(conservation_array[dr_zone_end + 1] >= CRASS_DEF_ZONE_EXT_CONS_CUT_OFF)
        {
            dr_zone_end++;
        }
        else
        {
        	break;
        }
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
        if(conservation_array[i] >= CRASS_DEF_COLLAPSED_CONS_CUT_OFF)
        {
            refined_DR_ends[i] = true;
            true_DR += consensus_array[i];
        }
        else
        {
            // possible collapsed cluster
            refined_DR_ends[i] = false;
#ifdef DEBUG
            logInfo("-------------", 5); 
            logInfo("Possible collapsed cluster at position: " << collapsed_pos << " (" << (dr_zone_start + collapsed_pos) << " || " << conservation_array[i] << ")", 5);
            logInfo("Base:  Count:  Cov:",5);
#endif
            float total_count = coverage_array[0][i] + coverage_array[1][i] + coverage_array[2][i] + coverage_array[3][i];
            
            for(int k = 0; k < 4; k++)
            {
#ifdef DEBUG
                logInfo("  " << alphabet[k] << "     " << coverage_array[k][i] << "      " << ((float)coverage_array[k][i]/total_count), 5);
#endif
                // check to make sure that each base is represented enough times
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
#ifdef DEBUG
                logInfo("   ...ignoring (FA)", 5);
#endif
                true_DR += consensus_array[i];
                refined_DR_ends[i] = true;
            }
            else
            {
                // is this seen at the DR level?
                std::map<char, int> collapsed_options2;
                DR_ClusterIterator dr_iter2 = clustered_DRs->begin();
                while (dr_iter2 != clustered_DRs->end()) 
                {
                    std::string tmp_DR = mStringCheck.getString(*dr_iter2);
                    if(-1 != DR_offset_map[*dr_iter2])
                    {
                        // check if the deciding character is within range of this DR
                    	// collapsed_pos + dr_zone_start gives the index in the ARRAY of the collapsed char
                    	// DR_offset_map[*dr_iter2] gives the start of the DR in the array
                    	// We need to check that collapsed_pos + dr_zone_start >= DR_offset_map[*dr_iter2] AND
                    	// that collapsed_pos < dr_zone_start - DR_offset_map[*dr_iter2] + tmp_DR.length()
                        //if(DR_offset_map[*dr_iter2] <= dr_zone_start && dr_zone_start < (DR_offset_map[*dr_iter2] + (int)(tmp_DR.length())) && collapsed_pos < (int)(tmp_DR.length()))
                    	if((collapsed_pos + dr_zone_start >= DR_offset_map[*dr_iter2]) && (collapsed_pos + dr_zone_start - DR_offset_map[*dr_iter2] < ((int)tmp_DR.length())))
                        {
                            // this is easy, we can compare based on this char only
                            char decision_char = tmp_DR[dr_zone_start - DR_offset_map[*dr_iter2] + collapsed_pos];
                            collapsed_options2[decision_char] = collapsed_options[decision_char];
                    		//logInfo("Index: " << (dr_zone_start - DR_offset_map[*dr_iter2] + collapsed_pos) << " dzs: " << dr_zone_start << " do: " << DR_offset_map[*dr_iter2] << " cp: " << collapsed_pos << " dl: " << ((int)tmp_DR.length()), 1);
                            //logInfo("Adding : " << decision_char << " to CO2", 1);
                        }
                    }
                    else
                    {
                    	logWarn("No offset for DR: " << tmp_DR, 1);
                    }
                    dr_iter2++;
                }
                
                if(2 > collapsed_options2.size())
                {
                    // in the case that the DR is collapsing at the very end of the zone,
                    // it may be because the spacers ahve a weird distribution of starting
                    // bases. We need to check this out here...
                    if(collapsed_pos == 0)
                    {
#ifdef DEBUG
                    	logInfo("   ...ignoring (RLO SS)", 5);
#endif
                    }
                    else if(collapsed_pos + dr_zone_start == dr_zone_end)
                    {
#ifdef DEBUG
                    	logInfo("   ...ignoring (RLO EE)", 5);
#endif
                    }
                    else
                    {
#ifdef DEBUG
                    	logInfo("   ...ignoring (RLO KK)", 5);
#endif
                        true_DR += consensus_array[i];
                        refined_DR_ends[i] = true;
                    }
                    collapsed_options.clear();
                }
                else
                {
                    // If it aint in collapsed_options2 it aint super fantastic!
                    collapsed_options.clear();
                    collapsed_options = collapsed_options2;
                    
                    // make the collapsed pos array specific and exit this loop
                    collapsed_pos += dr_zone_start;
                    i = dr_zone_end + 1;
                }
            }
        }
    }
    
    // check to make sure that the DR is not just some random long RE
    if((unsigned int)(true_DR.length()) > mOpts->highDRsize)
    {
        // probably a dud. throw it out
        // free the memory and clean up
        if(NULL != clustered_DRs)
        {
        	delete clustered_DRs;
        	clustered_DRs = NULL;
        	mDR2GIDMap[GID] = NULL;
        }

        logInfo("Killed: {" << true_DR << "} cause' it was too long", 1);

        //++++++++++++++++++++++++++++++++++++++++++++++++
        // clean up the mess we made
        
        if(NULL != consensus_array)
        	delete[] consensus_array;
        if(NULL != conservation_array)
        	delete[] conservation_array;
        if(coverage_array != NULL)
        {
    		for(int i = 0; i < 4; i++)
    		{
    			if(NULL != coverage_array[i])
    				delete[] coverage_array[i];
    		}
    		delete[] coverage_array;
    		coverage_array = NULL;
        }
        return false;
    }
    if(((unsigned int)(true_DR.length()) < mOpts->lowDRsize) && (collapsed_options.size() == 0))
    {
        // probably a dud. throw it out
        // free the memory and clean up
        if(NULL != clustered_DRs)
        {
        	delete clustered_DRs;
        	clustered_DRs = NULL;
        	mDR2GIDMap[GID] = NULL;
        }
        
        logInfo("Killed: {" << true_DR << "} cause' the consensus was too short...", 1);

        //++++++++++++++++++++++++++++++++++++++++++++++++
        // clean up the mess we made
        
        if(NULL != consensus_array)
        	delete[] consensus_array;
        if(NULL != conservation_array)
        	delete[] conservation_array;
        if(coverage_array != NULL)
        {
    		for(int i = 0; i < 4; i++)
    		{
    			if(NULL != coverage_array[i])
    				delete[] coverage_array[i];
    		}
    		delete[] coverage_array;
    		coverage_array = NULL;
        }
        return false;
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // print out the consensus array 
    
    if(0 == collapsed_options.size())
    {
        // update the DR start and ends
        int diffs = dr_zone_end - dr_zone_start + 1 - (int)true_DR.length();
        while(0 < diffs)
        {
            // we need to update the start or end
            if(!refined_DR_ends[dr_zone_end])
            {
                dr_zone_end--;
                diffs--;
            }
            if(0 < diffs)
            {
				if(!refined_DR_ends[dr_zone_start])
				{
					dr_zone_start++;
					diffs--;
				}
            }
        }
        
        if(isLogging(3))
        {
        	int show_xtra = 4;
        	int print_start = dr_zone_start - show_xtra;
        	int print_end = dr_zone_end + show_xtra;
        	stringstream ss;
            ss << std::endl << "%, ";
            for(int i = print_start; i <= print_end; i++)
            {
            	if(i == dr_zone_start)
            		ss << "|,"; 
                ss << conservation_array[i] << ", ";
            	if(i == dr_zone_end)
            		ss << "|,"; 
            }
            
            for(int j = 0; j < 4; j++)
            {
                ss << std::endl;
                ss << alphabet[j] << ", ";
                for(int i = print_start; i <= print_end; i++)
                {
                	if(i == dr_zone_start)
                		ss << "|,"; 
                    ss << coverage_array[j][i] << ", ";
                	if(i == dr_zone_end)
                		ss << "|,"; 
                }
            } 
            ss << std::endl << "$, ";
            for(int i = print_start; i <= print_end; i++)
            {
            	if(i == dr_zone_start)
            		ss << "|,"; 
                ss << consensus_array[i] << ", ";
            	if(i == dr_zone_end)
            		ss << "|,"; 
            }
            logInfo(ss.str(), 3);
        }
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // clean up the mess we made
    
    if(NULL != consensus_array)
    {
        delete[] consensus_array;
    }
    if(NULL != conservation_array)
    {
        delete[] conservation_array;
    }
    if(coverage_array != NULL)
    {
		for(int i = 0; i < 4; i++)
		{
			if(NULL != coverage_array[i])
            {
                delete[] coverage_array[i];
            }
		}
		delete[] coverage_array;
		coverage_array = NULL;
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
            mDR2GIDMap[group] = new DR_Cluster;
            coll_char_to_GID_map[co_iter->first] = group;
            logInfo("Mapping \""<< co_iter->first <<"\" to group: " << group, 1);
            co_iter++;
        }
        
        dr_iter = clustered_DRs->begin();
        while (dr_iter != clustered_DRs->end()) 
        {
            std::string tmp_DR = mStringCheck.getString(*dr_iter);
            if(-1 != DR_offset_map[*dr_iter])
            {
                // check if the deciding character is within range of this DR
                if(DR_offset_map[*dr_iter] <= collapsed_pos && collapsed_pos < (DR_offset_map[*dr_iter] + (int)(tmp_DR.length())))
                {
                    // this is easy, we can compare based on this char only
                    char decision_char = tmp_DR[collapsed_pos - DR_offset_map[*dr_iter]];
                    (mDR2GIDMap[ coll_char_to_GID_map[ decision_char ] ])->push_back(*dr_iter);
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
                        StartStopListIterator ss_iter = (*read_iter)->begin();
                        while(ss_iter != (*read_iter)->end())
                        {
                            int within_read_dec_pos = *ss_iter + dec_diff;
                            if(within_read_dec_pos > 0 && within_read_dec_pos < (int)(*read_iter)->getSeqLength())
                            {
                                char decision_char = (*read_iter)->getSeqCharAt(within_read_dec_pos);
                                
                                // it must be one of the collapsed options!
                                if(collapsed_options.find(decision_char) != collapsed_options.end())
                                {
                                	forms_map[decision_char] = NULL;
                                	break;
                                }
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
                                StartStopListIterator ss_iter = (*read_iter)->begin();
                                while(ss_iter != (*read_iter)->end())
                                {
                                    int within_read_dec_pos = *ss_iter + dec_diff;
                                    if(within_read_dec_pos > 0 && within_read_dec_pos < (int)(*read_iter)->getSeqLength())
                                    {
                                        char decision_char = (*read_iter)->getSeqCharAt(within_read_dec_pos);
                                        // it must be one of the collapsed options!
                                        if(forms_map.find(decision_char) != forms_map.end())
                                        {
                                        	(mDR2GIDMap[ coll_char_to_GID_map[ decision_char ] ])->push_back(*dr_iter);
                                            break_out = true;
                                            break;
                                        }
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
                            // Something is wrong!
#ifdef DEBUG                        	
                            logWarn("No reads fit the form: " << tmp_DR, 8);
#endif
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
                                StringToken st = mStringCheck.addString(tmp_DR);
                                mReads[st] = new ReadList();
                                // make sure we know which readlist is which
                                forms_map[fm_iter->first] = mReads[st];
                                // put the new dr_token into the right cluster
                                (mDR2GIDMap[ coll_char_to_GID_map[ fm_iter->first ] ])->push_back(st);
                                
                                // next!
                                fm_iter++;
                            }
                            
                            // put the correct reads on the correct readlist
                            read_iter = mReads[*dr_iter]->begin();
                            while (read_iter != mReads[*dr_iter]->end()) 
                            {
                                StartStopListIterator ss_iter = (*read_iter)->begin();
                                while(ss_iter != (*read_iter)->end())
                                {
                                    int within_read_dec_pos = *ss_iter + dec_diff;
                                    if(within_read_dec_pos > 0 && within_read_dec_pos < (int)(*read_iter)->getSeqLength())
                                    {
                                        char decision_char = (*read_iter)->getSeqCharAt(within_read_dec_pos);
                                        
                                        // needs to be a form we've seen before!
                                        if(forms_map.find(decision_char) != forms_map.end())
                                        {
											// push this readholder onto the correct list
											(forms_map[decision_char])->push_back(*read_iter);
											
											// make the original pointer point to NULL so we don't delete twice
											*read_iter = NULL;
											
											break;
                                        }
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
        if(NULL != clustered_DRs)
        {
        	delete clustered_DRs;
        	clustered_DRs = NULL;
        	mDR2GIDMap[GID] = NULL;
        }
        
        logInfo("Calling the parser recursively", 4);
        
        // call this baby recursively with the new clusters
        std::map<char, int>::iterator cc_iter = coll_char_to_GID_map.begin();
        while(cc_iter != coll_char_to_GID_map.end())
        {
            parseGroupedDRs(cc_iter->second, nTopKmers, mDR2GIDMap[cc_iter->second], nextFreeGID);
            cc_iter++;
        }
    }
    else
    {
        //++++++++++++++++++++++++++++++++++++++++++++++++
        // repair all the startstops for each read in this group
        //
        // This function is recursive, so we'll only get here when we have found exactly one DR
        
        logInfo("Found DR: " << true_DR, 2);
        mTrueDRs[GID] = true_DR;
		DR_ClusterIterator drc_iter = (mDR2GIDMap[GID])->begin();
		while(drc_iter != (mDR2GIDMap[GID])->end())
		{
			if (DR_offset_map[*drc_iter] == -1 ) 
            {
                logError("Repeat "<< *drc_iter<<" in Group "<<GID <<" has no offset in DR_offset_map");
            } 
            else 
            {
                // go through each read
                ReadListIterator read_iter = mReads[*drc_iter]->begin();
                while (read_iter != mReads[*drc_iter]->end()) 
                {
                    (*read_iter)->updateStartStops((DR_offset_map[*drc_iter] - dr_zone_start), &true_DR, mOpts);
                    read_iter++;
                }
            }

			drc_iter++;
		}
    }
    
    return true;
}
int WorkHorse::numberOfReadsInGroup(DR_Cluster * currentGroup)
{
    DR_ClusterIterator grouped_drs_iter = currentGroup->begin();
    size_t number_of_reads_in_group = 0;
    while (grouped_drs_iter != currentGroup->end()) 
    {
        number_of_reads_in_group += mReads[*grouped_drs_iter]->size();
        ++grouped_drs_iter;
    }
    return (int)number_of_reads_in_group;
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

bool WorkHorse::getNMostAbundantKmers(std::vector<std::string>& mostAbundantKmers, int num2Get, std::map<std::string, int> * kmer_CountMap)
{
    std::string top_kmer;    
    std::map<std::string, bool> top_kmer_map;
    
    if ((int)(kmer_CountMap->size()) < num2Get) 
    {
        return false;
    } 
    else 
    {
        for (int i = 1; i <= num2Get; i++) 
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
        std::map<std::string, bool>::iterator tkm_iter = top_kmer_map.begin();
        while(tkm_iter != top_kmer_map.end())
        {
            //std::cout<<tkm_iter->first<<std::endl;
            mostAbundantKmers.push_back(tkm_iter->first);
            tkm_iter++;
        }
        return true;
    }
}

bool WorkHorse::clusterDRReads(StringToken DRToken, int * nextFreeGID, std::map<std::string, int> * k2GIDMap, std::map<int, std::map<std::string, int> * > * groupKmerCountsMap)
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
        kmers[i][CRASS_DEF_KMER_SIZE] = '\0';
        
        std::string tmp_str(kmers[i]);
        tmp_str = laurenize(tmp_str);
        
        // see if this guy has been counted!
        if(local_kmer_CountMap.find(tmp_str) == local_kmer_CountMap.end())
        {
            local_kmer_CountMap[tmp_str] = 1;
        }
        else
        {
            local_kmer_CountMap[tmp_str]++;
        }
        
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
        mGroupMap[group] = true;
        mDR2GIDMap[group] = new DR_Cluster;
        
        // we need a new kmer counter for this group
        (*groupKmerCountsMap)[group] = new std::map<std::string, int>;
    }
    
    // we need to record the group for this mofo!
    mDR2GIDMap[group]->push_back(DRToken);
    
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
    
    return true;
    
}

//**************************************
// spacer graphs
//**************************************
int WorkHorse::makeSpacerGraphs(void)
{
	//-----
	// build the spacer graphs
	//
    // go through the DR2GID_map and make all reads in each group into nodes
    DR_ListIterator dr_iter = mDRs.begin();
    while(dr_iter != mDRs.end())
    {
        if(NULL != dr_iter->second)
        {
        	logInfo("Making spacer graph for DR: " << dr_iter->first, 1);
        	if((dr_iter->second)->buildSpacerGraph())
        		return 1;
        }
        dr_iter++;
    }
    return 0;
}

int WorkHorse::cleanSpacerGraphs(void)
{
	//-----
	// clean the spacer graphs
	//
    // go through the DR2GID_map and make all reads in each group into nodes
    DR_ListIterator dr_iter = mDRs.begin();
    while(dr_iter != mDRs.end())
    {
        if(NULL != dr_iter->second)
        {
        	logInfo("Cleaning spacer graph for DR: " << dr_iter->first, 1);
        	//(dr_iter->second)->printAllSpacers();
        	if((dr_iter->second)->cleanSpacerGraph())
        		return 1;
        }
        dr_iter++;
    }
    return 0;
}

int WorkHorse::generateFlankers(void)
{
	//-----
	// Wrapper for graph cleaning
	//
	// create a spacer dictionary
	logInfo("Detecting Flanker sequences", 1);
	DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(NULL != drg_iter->second)
		{            
			if (NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                (mDRs[mTrueDRs[drg_iter->first]])->generateFlankers();
		    }
        }
		drg_iter++;
	}
	return 0;
}
//**************************************
// contig making
//**************************************
int WorkHorse::splitIntoContigs(void)
{
	//-----
	// split all groups into contigs
	//
    // go through the DR2GID_map and make all reads in each group into nodes
    DR_ListIterator dr_iter = mDRs.begin();
    while(dr_iter != mDRs.end())
    {
        if(NULL != dr_iter->second)
        {
        	logInfo("Making spacer contigs for DR: " << dr_iter->first, 1);
        	if((dr_iter->second)->splitIntoContigs())
        		return 1;
        }
        dr_iter++;
    }
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
    std::string kmer_lookup_name = mOutFileDir + CRASS_DEF_DEF_KMER_LOOKUP_EXT;
    std::string patterns_lookup_name = mOutFileDir + CRASS_DEF_DEF_PATTERN_LOOKUP_EXT;
    std::string spacer_lookup_name = mOutFileDir + CRASS_DEF_DEF_SPACER_LOOKUP_EXT;
    
    // Write!
    writeLookupToFile(kmer_lookup_name, kmerLookup);  
    writeLookupToFile(patterns_lookup_name, patternsLookup);
    writeLookupToFile(spacer_lookup_name, spacerLookup);
}


int WorkHorse::dumpReads(DR_Cluster_Map * DR2GID_map, bool split)
{
    //-----
    // Print the reads from one cluster to a file...
    //
	logInfo("Dumping reads", 1);
    DR_Cluster_MapIterator drg_iter = DR2GID_map->begin();
    while (drg_iter != DR2GID_map->end()) 
    {
        // make sure that our cluster is real
        if (drg_iter->second != NULL) 
        {
            if(NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                (mDRs[mTrueDRs[drg_iter->first]])->dumpReads((mOpts->output_fastq +  "Group_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first] + ".fa").c_str(), false, split);
            }
        }
        drg_iter++;
    }
	return 0;
}

int WorkHorse::dumpSpacers(void)
{
	//-----
	// Wrapper for graph cleaning
	//
	// create a spacer dictionary
	logInfo("Dumping spacers", 1);
	DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(NULL != drg_iter->second)
		{            
			if (NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                (mDRs[mTrueDRs[drg_iter->first]])->dumpSpacerDict(mOpts->output_fastq + "Group_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first] + ".spacers", false);
		    }
        }
		drg_iter++;
	}
	return 0;
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

int WorkHorse::renderDebugGraphs(void)
{
	//-----
	// Print the debug graph
	//
    // use the default name
	return renderDebugGraphs("Group_");
}

int WorkHorse::renderDebugGraphs(std::string namePrefix)
{
	//-----
	// Print the debug graph
	//
	// go through the DR2GID_map and make all reads in each group into nodes
    logInfo("Rendering debug graphs" , 1);
    
    DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
    while(drg_iter != mDR2GIDMap.end())
    {
        if(NULL != drg_iter->second)
        {            
            if (NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                std::ofstream graph_file;
                std::string graph_file_prefix = mOpts->output_fastq + namePrefix + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first];
                std::string graph_file_name = graph_file_prefix + "_debug.gv";
                graph_file.open(graph_file_name.c_str());
                if (graph_file.good()) 
                {
                    mDRs[mTrueDRs[drg_iter->first]]->printDebugGraph(graph_file, mTrueDRs[drg_iter->first], false, false, false);
#if RENDERING
                    if (!mOpts->noRendering) 
                    {
                        // create a command string and call neato to make the image file
                        std::string cmd = "neato -Teps " + graph_file_name + " > "+ graph_file_prefix + ".eps";
                        if (system(cmd.c_str()))
                        {
                            logError("Problem running neato when rendering debug graphs");
                        }
                    }
#endif
                } 
                else 
                {
                    logError("Unable to create graph output file "<<graph_file_name);
                }
                graph_file.close();
            }
        }
        drg_iter++;
    }
    return 0;
}

int WorkHorse::renderSpacerGraphs(void)
{
	//-----
	// Print the cleaned? spacer graph
	//
    // use the default name
	return renderSpacerGraphs("Spacers_");
}

int WorkHorse::renderSpacerGraphs(std::string namePrefix)
{
	//-----
	// Print the cleaned? spacer graph
	//
	// go through the DR2GID_map and make all reads in each group into nodes
    logInfo("Rendering spacer graphs" , 1);
    
    // make a single file with all of the keys for the groups
    std::ofstream key_file;
    
    std::stringstream key_file_name;
    key_file_name << mOpts->output_fastq<<PACKAGE_NAME << "_"<<mTimeStamp<<"_keys.gv";
    key_file.open(key_file_name.str().c_str());

    if (!key_file) 
    {
        logError("Cannot open the key file");
        return 1;
    }

    gvGraphHeader(key_file, "Keys");
    DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
    while(drg_iter != mDR2GIDMap.end())
    {
        if(NULL != drg_iter->second)
        {            
            if(NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                std::ofstream graph_file;

                
                std::string graph_file_prefix = mOpts->output_fastq + namePrefix + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first];
                std::string graph_file_name = graph_file_prefix + "_spacers.gv";
                graph_file.open(graph_file_name.c_str());
                if (graph_file.good()) 
                {
                    mDRs[mTrueDRs[drg_iter->first]]->printSpacerGraph(graph_file, mTrueDRs[drg_iter->first], mOpts->longDescription, mOpts->showSingles);
                    mDRs[mTrueDRs[drg_iter->first]]->printSpacerKey(key_file, 10, namePrefix + to_string(drg_iter->first));
#if RENDERING
                    if (!mOpts->noRendering) 
                    {
                        // create a command string and call graphviz to make the image file
                        std::string cmd = mOpts->layoutAlgorithm + " -Teps " + graph_file_name + " > "+ graph_file_prefix + ".eps";
                        if(system(cmd.c_str()))
                        {
                            logError("Problem running "<<mOpts->layoutAlgorithm<<" when rendering spacer graphs");
                            return 1;
                        }
                    }
#endif
                } 
                else 
                {
                    logError("Unable to create graph output file "<<graph_file_name);
                    return 1;
                }
                graph_file.close();
            }
        }
        drg_iter++;
    }
    gvGraphFooter(key_file);
    key_file.close();
    return 0;
}

bool WorkHorse::checkFileOrError(const char * fileName)
{
    try {
        // Test to see if the file is ok.
        struct stat inputDirStatus;
        int xStat = stat(fileName, &inputDirStatus);
        // stat failed
        switch (xStat) 
        {
            case -1:
            {
                switch (errno)
                {
                    case ENOENT:
                    {
                        throw ( std::runtime_error("Path to file does not exist, or path is an empty string.") );
                        break;
                    }
                    case ELOOP:
                    {
                        throw ( std::runtime_error("Too many symbolic links encountered while traversing the path to file."));
                        break;
                    }
                    case EACCES:
                    {
                        throw ( std::runtime_error("You do not have permission to access the file."));
                        break;
                    }
                    default:
                    {
                        throw (std::runtime_error("An error occured when reading the file"));
                        break;
                    }
                }
                break;
            }
            default:
            {
                return true;
                break;
            }
        }
    } catch (std::exception& e) {
        std::cerr << e.what()<<std::endl;
        logError(e.what());
        return false;

    }
}


bool WorkHorse::printXML(std::string namePrefix)
{
	// print all the assembly gossip to XML
	namePrefix += CRASS_DEF_CRISPR_EXT;
	logInfo("Writing XML output to \"" << namePrefix << "\"", 1);
	
    CrassXML * xml_doc = new CrassXML();
    int error_num;
    xercesc::DOMElement * root_element = xml_doc->createDOMDocument(CRASS_DEF_ROOT_ELEMENT, CRASS_DEF_XML_VERSION, error_num);
    
    if (root_element && !error_num) 
    {
        // go through the node managers and print the group info 
        // print all the inside information
        DR_Cluster_MapIterator drg_iter =  mDR2GIDMap.begin();
        while (drg_iter != mDR2GIDMap.end()) 
        {
            // make sure that our cluster is real
            if (drg_iter->second != NULL) 
            {
                if(NULL != mDRs[mTrueDRs[drg_iter->first]])
                {
                    std::string gid_as_string = "G" + to_string(drg_iter->first);
                    xercesc::DOMElement * group_elem = xml_doc->addGroup(gid_as_string, mTrueDRs[drg_iter->first], root_element);
                    
                    
                    /*
                     * <data> section
                     */
                    addDataToDOM(xml_doc, group_elem, drg_iter->first);

                    /*
                     * <metadata> section
                     */
                    addMetadataToDOM(xml_doc, group_elem, drg_iter->first);
                    
                    /*
                     * <assembly> section
                     */
                    xercesc::DOMElement * assem_elem = xml_doc->addAssembly(group_elem);
                    (mDRs[mTrueDRs[drg_iter->first]])->printAssemblyToDOM(xml_doc, assem_elem, false);
                    
                }
            }
            drg_iter++;
        }
        xml_doc->printDOMToFile(namePrefix);
    } 
    else
    {
        return 1;
    }
    delete xml_doc;
	return 0;
}

bool WorkHorse::addDataToDOM(CrassXML * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    try 
    {
        xercesc::DOMElement * data_elem = xmlDoc->addData(groupElement);
        if ((mDRs[mTrueDRs[groupNumber]])->haveAnyFlankers()) {
            std::cout<<"got some flankers"<<std::endl;
            xmlDoc->createFlankers(data_elem);
        }
        
        for (xercesc::DOMElement * currentElement = data_elem->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) 
        {
            if( xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->getDrs()))
            {
                // TODO: current implementation in Crass only supports a single DR for a group
                // in the future this will change, but for now ok to keep as a constant
                std::string drid = "DR1";
                xmlDoc->addDirectRepeat(drid, mTrueDRs[groupNumber], currentElement);
            }
            else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->getSpacers()))
            {
                // print out all the spacers for this group
                (mDRs[mTrueDRs[groupNumber]])->addSpacersToDOM(xmlDoc, currentElement, false);
                
            }
            else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->getFlankers()))
            {
                // should only get in here if there are flankers for the group

                // print out all the flankers for this group
                (mDRs[mTrueDRs[groupNumber]])->addFlankersToDOM(xmlDoc, currentElement, false);

            }
        }
    }
    catch( xercesc::XMLException& e )
    {
        char* message = xercesc::XMLString::transcode( e.getMessage() );
        std::ostringstream errBuf;
        errBuf << "Error parsing file: " << message << std::flush;
        xercesc::XMLString::release( &message );
        return 1;
    }
    return 0;
}

bool WorkHorse::addMetadataToDOM(CrassXML * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    std::stringstream notes;
    notes << PACKAGE_NAME <<" ("<<PACKAGE_VERSION<<") run on "<<mTimeStamp<<" with command: ";
    notes <<mCommandLine;
    xercesc::DOMElement * metadata_elem = xmlDoc->addMetaData(notes.str(), groupElement);
    
    std::string file_name;
    // add in files if they exist
    if (!mOpts->logToScreen) 
    {
        // we whould have a log file
        file_name = mOpts->output_fastq + PACKAGE_NAME + "." + mTimeStamp + ".log";
        if (checkFileOrError(file_name.c_str())) 
        {
            xmlDoc->addFileToMetadata("log", file_name, metadata_elem);
        }
        else
        {
            logError("Could not find the log file at "<<file_name<<" but I think it should be there... wierd");
        }
    }
    
    
#ifdef DEBUG
    // check for debuging .gv files
    file_name = mOpts->output_fastq + "Group_"; 
    std::string file_sufix = to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + "_debug.gv";
    if (checkFileOrError((file_name + file_sufix).c_str())) 
    {
        xmlDoc->addFileToMetadata("data", (file_name + file_sufix), metadata_elem);
    } 
    else 
    {
        logError("Could not find the Debug .gv file at "<< file_name << file_sufix <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }
    
    // and now for the cleaned .gv
    file_name = mOpts->output_fastq + "Clean_";
    if (checkFileOrError((file_name + file_sufix).c_str())) 
    {
        xmlDoc->addFileToMetadata("data", (file_name + file_sufix), metadata_elem);
    } 
    else 
    {
        logError("Could not find the Debug .gv file at "<<file_name << file_sufix <<" for group " <<groupNumber<<", but I think it should be there... wierd");

    }
    
#endif
    
#ifdef RENDERING
    // check for image files
#ifdef DEBUG
    file_name = mOpts->output_fastq + "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
    if (checkFileOrError(file_name.c_str())) 
    {
        xmlDoc->addFileToMetadata("image", file_name, metadata_elem);
    } 
    else 
    {
        logError("Could not find the Debug .eps file at "<<file_name <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }
    
    file_name = mOpts->output_fastq + "Clean_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
    
    if (checkFileOrError(file_name.c_str())) 
    {
        xmlDoc->addFileToMetadata("image", file_name, metadata_elem);
    } 
    else 
    {
        logError("Could not find the cleaned debug .eps file at "<<file_name <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }
#endif // DEBUG
    
    file_name = mOpts->output_fastq + "Spacers_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
    if (checkFileOrError(file_name.c_str())) 
    {
        xmlDoc->addFileToMetadata("image", file_name, metadata_elem);
    } 
    else 
    {
        logError("Could not find the Spacer .eps file at "<<file_name <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }

#endif
    
    // check the sequence file
    file_name = mOpts->output_fastq +  "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".fa";
    if (checkFileOrError(file_name.c_str())) 
    {
        xmlDoc->addFileToMetadata("sequence", file_name, metadata_elem);
    } 
    else 
    {
        logError("Could not find the fasta file at "<<file_name <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }
    
    
    // check the spacer dictionary 
    file_name = mOpts->output_fastq +  "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".spacers";
    if (checkFileOrError(file_name.c_str())) 
    {
        xmlDoc->addFileToMetadata("data", file_name, metadata_elem);
    } 
    else 
    {
        logError("Could not find the spacer dictionary file at "<<file_name <<" for group " <<groupNumber<<", but I think it should be there... wierd");
    }
    return 0;
    
}



