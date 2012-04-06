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
#include <unistd.h>

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
#include "Exception.h"

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
    // Not needed any more since we have the crispr file
//	if(dumpSpacers())
//	{
//        logError("FATAL ERROR: dumpSpacers failed");
//        return 8;
//	}
    
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
        try {
            mAveReadLength += decideWhichSearch(seq_iter->c_str(), *mOpts, &mReads, &mStringCheck, patterns_lookup, reads_found);
            logInfo("Finished file: " << *seq_iter, 1);

        } catch (crispr::exception& e) {
            std::cerr<<e.what()<<std::endl;
            return 1;
        }        
        seq_iter++;
    }
    mAveReadLength /= (double)seqFiles.size();
    if (patterns_lookup.size() > 0) 
    {
        logInfo("Begining Second iteration through files to recruit singletons", 2);
        std::cout<<"["<<PACKAGE_NAME<<"_patternFinder]: " << patterns_lookup.size() << " patterns."<<std::endl;
        seq_iter = seqFiles.begin();
        while (seq_iter != seqFiles.end()) {
            
            logInfo("Parsing file: " << *seq_iter, 1);
            
            try {
                findSingletons(seq_iter->c_str(), *mOpts, patterns_lookup, reads_found, &mReads, &mStringCheck);
            } catch (crispr::exception& e) {
                std::cerr<<e.what()<<std::endl;
                return 1;
            }
            seq_iter++;
        }
    }
    logInfo("Searching complete. " << mReads.size()<<" direct repeat variants have been found", 1);
    logInfo("Number of reads found so far: "<<this->numOfReads(), 2);
    
    // There will be an abundance of forms for each direct repeat.
    // We needs to do somes clustering! Then trim and concatenate the direct repeats
    try {
        if (mungeDRs())
        {
            logError("Wierd stuff happend when trying to get the 'true' direct repeat");            
            return 1;
        }
    } catch(crispr::exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
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
    std::cout<<'['<<PACKAGE_NAME<<"_graphBuilder]: "<<mTrueDRs.size()<<" putative CRISPRs found!"<<std::endl;
    //MI std::cout<<'['<<PACKAGE_NAME<<"_graphBuilder]: "<<std::flush;
    while(drg_iter != mDR2GIDMap.end())
    {
        if(NULL != drg_iter->second)
        {            
#ifdef DEBUG
            logInfo("Creating NodeManager "<<drg_iter->first, 6);
#endif
            //MI std::cout<<'['<<drg_iter->first<<','<<mTrueDRs[drg_iter->first]<<std::flush;
            mDRs[mTrueDRs[drg_iter->first]] = new NodeManager(mTrueDRs[drg_iter->first], mOpts);
            //MI std::cout<<'.'<<std::flush;
            DR_ClusterIterator drc_iter = (drg_iter->second)->begin();
            while(drc_iter != (drg_iter->second)->end())
            {
                // go through each read
            	//MI std::cout<<'|'<<std::flush;
                ReadListIterator read_iter = mReads[*drc_iter]->begin();
                while (read_iter != mReads[*drc_iter]->end()) 
                {
                	//MI std::cout<<'.'<<std::flush;
                    mDRs[mTrueDRs[drg_iter->first]]->addReadHolder(*read_iter);
                    read_iter++;
                }
                drc_iter++;
            }
            //MI std::cout<<"],"<<std::flush;
        }
        drg_iter++;
    }
    //MI std::cout<<std::endl;
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
            if (NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                if((mDRs[mTrueDRs[drg_iter->first]])->getSpacerCount(false) < mOpts->covCutoff) 
                {
                    logInfo("Deleting NodeManager "<<drg_iter->first<<" as it contained less than "<<mOpts->covCutoff<<" attached spacers",4);
                    delete mDRs[mTrueDRs[drg_iter->first]];
                    mDRs[mTrueDRs[drg_iter->first]] = NULL;
                }
            }
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
    std::cout<<'['<<PACKAGE_NAME<<"_clusterCore]: "<<mReads.size()<<" variants mapped to "<<mDR2GIDMap.size()<<" clusters"<<std::endl;

    if(isLogging(4))
    {
        DR_Cluster_MapIterator dcg_iter = mDR2GIDMap.begin();
        while(dcg_iter != mDR2GIDMap.end())
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
            dcg_iter++;
        } 
    }

    logInfo("Reducing list of potential DRs (2): Cluster refinement and true DR finding", 1);
    
    // go through all the counts for each group
    std::map<int, std::map<std::string, int> * >::iterator group_count_iter = group_kmer_counts_map.begin();
    while(group_count_iter != group_kmer_counts_map.end())
    {
        if(NULL != mDR2GIDMap[group_count_iter->first])
        {
            // it's real, so parse this group
            // get the N top kmers
            std::vector<std::string> n_top_kmers;
            int num_mers_found = getNMostAbundantKmers((int)(mDR2GIDMap[group_count_iter->first])->size(), n_top_kmers, CRASS_DEF_NUM_KMERS_4_MODE, group_count_iter->second);
            if (0 != num_mers_found) 
            {
                // a return value of false indicates that this function has deleted clustered_DRs
            	//std::cout << "found " << num_mers_found << " : " << CRASS_DEF_NUM_KMERS_4_MODE << std::endl; 
                parseGroupedDRs(num_mers_found, group_count_iter->first, &n_top_kmers, &next_free_GID);
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

bool WorkHorse::findMasterDR(int GID, std::vector<std::string> * nTopKmers, StringToken * masterDRToken, std::string * masterDRSequence)
{
	//-----
	// Identify a master DR
	// 
	// Updates the values in masterDRToken and  masterDRSequence and returns true
	// otherwise deletes the memory pointed at by clustered_DRs and returns false
	//
	//
    // to store our DR which has all XX kmers
	logInfo("Identifying a master DR", 1);
    int master_read_count = 0;
    
    // these are needed for the call to is kmer present but we don't actually need to values!
    bool disp_rc;
    int disp_pos;
    
    // go through the DRs in this cluster, we'd like to find one which has all kmers in it...
    // moreover, we need to get the one with all XX and the most reads
    DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
    while (dr_iter != (mDR2GIDMap[GID])->end()) 
    {
        std::string tmp_DR = mStringCheck.getString(*dr_iter);
        bool got_all_mode_mers = true;
        std::vector<std::string>::iterator n_top_iter = (*nTopKmers).begin();
        while(n_top_iter != (*nTopKmers).end())
        {
            if(!isKmerPresent(&disp_rc, &disp_pos, *n_top_iter, &tmp_DR))
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
                *masterDRSequence = mStringCheck.getString(*dr_iter);
                *masterDRToken = *dr_iter;
                logInfo("Identified: " << *masterDRSequence << " (" << *masterDRToken << ") as a master potential DR", 4);
                return true;
            }
        }
        
        // otherwise keep searching
        dr_iter++;
    }
    
    if(*masterDRToken == -1)
    {
        // probably a dud. throw it out
        // free the memory and clean up
        logInfo("Could not identify a master DR", 4);
        if(NULL != mDR2GIDMap[GID])
        {
        	delete mDR2GIDMap[GID];
        	mDR2GIDMap[GID] = NULL;
            mDR2GIDMap.erase(GID);
        }
    }
    
    return false;
}

bool WorkHorse::populateCoverageArray(int numMers4Mode, int GID, std::string master_DR_sequence, StringToken master_DR_token, std::map<StringToken, int> * DR_offset_map, int * dr_zone_start, int * dr_zone_end, std::vector<std::string> * nTopKmers, int ** coverage_array, int * kmer_positions_DR, bool * kmer_rcs_DR, int * kmer_positions_DR_master, bool * kmer_rcs_DR_master, int * kmer_positions_ARRAY)
{
	//-----
	// Use the data structures initialised in parseGroupedDRs
	// Load all the reads into the consensus array
	//
	logInfo("Populating consensus array", 1);

	bool first_run = true;          // we only need to do this once
	int array_len = ((int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength > CRASS_DEF_MIN_CONS_ARRAY_LEN) ? (int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength : CRASS_DEF_MIN_CONS_ARRAY_LEN;

	// chars we luv!
    char alphabet[4] = {'A', 'C', 'G', 'T'};

    // First we add the master DR into the arrays 
    ReadListIterator read_iter = mReads[master_DR_token]->begin();
    while (read_iter != mReads[master_DR_token]->end()) 
    {
        // don't care about partials
        int dr_start_index = 0;
        int dr_end_index = 1;

        // Find the DR which is the reported length. 
        // NOTE: If you -1 from the master_DR_sequence.length() this loop will throw an error as it will never be true
        while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) != ((int)(master_DR_sequence.length()) - 1))
        {
            dr_start_index += 2;
            dr_end_index += 2;
        }
        
        //  This if is to catch some weird-ass scenario, if you get a report that everything is wrong, then you've most likely
        // corrupted memory somewhere!
        if(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == ((int)(master_DR_sequence.length()) - 1))
        {
            // the start of the read is 
            int this_read_start_pos = kmer_positions_ARRAY[0] - (*read_iter)->startStopsAt(dr_start_index) - kmer_positions_DR[0] ;
            if(first_run)
            {
                *dr_zone_start =  this_read_start_pos + (*read_iter)->startStopsAt(dr_start_index);
                *dr_zone_end =  this_read_start_pos + (*read_iter)->startStopsAt(dr_end_index);
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
                    default:
                        index = 3;
                        break;
                }
				int index_b = i+this_read_start_pos; 
				if((index_b) >= array_len)
				{
					logError("The consensus/coverage arrays are too short. Consider changing CRASS_DEF_MIN_CONS_ARRAY_LEN to something larger and re-compiling");
				}
				if((index_b) < 0)
				{
					logError("AARRGGHH!!! " << i << " : " << this_read_start_pos << " : " << index_b << " : " << kmer_positions_ARRAY[0]  << " - " << (*read_iter)->startStopsAt(dr_start_index) << " - " << kmer_positions_DR[0]);
				}

				coverage_array[index][index_b]++;
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
    DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
    while (dr_iter != (mDR2GIDMap[GID])->end()) 
    {
        // get the string for this mofo
        std::string tmp_DR = mStringCheck.getString(*dr_iter);
        
        // we've already done the master DR
        if(master_DR_token != *dr_iter)
        {

            // set this guy to -1 for now
            (*DR_offset_map)[*dr_iter] = -1;
            
            // this is a DR we have yet to add to the coverage array
            // First we need to find the positions of the kmers in this DR

            for(int i = 0; i < numMers4Mode; i++)
            {
                isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), (*nTopKmers)[i], &tmp_DR);
                
            }            
            // we need to have at least half of the mode k-mers present to continue
            // Where they are not found, the position gets set to -1
            int num_neg1 = 0;
            for(int i = 0; i < numMers4Mode; i++)
            {
                if(-1 == (kmer_positions_DR)[i])
                    num_neg1++; 
            }
            int numMers4Mode_half = numMers4Mode / 2;
            //std::cout<<numMers4Mode_half<<" : " <<numMers4Mode<< " : "<<num_neg1<<std::endl;
            if(num_neg1 < numMers4Mode_half)
            {
                // all is good! ...so far
                // now we can determine if the DR is the reverse complement...
                int num_agree = 0;
                int num_disagree = 0;
                for(int i = 0; i < numMers4Mode; i++)
                {
                    if((kmer_rcs_DR_master)[i] == (kmer_rcs_DR)[i])
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
                        (*DR_offset_map)[*dr_iter] = -1;
                    }
                    
                    // TODO: Maybe we should use the smith-waterman here
                    // So now that all of the variants are in the same 
                    // orientation, perform a pair-wise alignment against the 
                    // master to get the offset in the array
                    //std::string master_dr_seq = mStringCheck.getString(master_DR_token);
                    //getOffsetAgainstMaster(master_dr_seq, tmp_DR);
                    
                    
                    // find the position of the kmers in the direct repeat
                    // store whether it was revcomp'd or not
                    for(int i = 0; i < numMers4Mode; i++)
                    {
                        isKmerPresent(((kmer_rcs_DR) + i), ((kmer_positions_DR) + i), (*nTopKmers)[i], &tmp_DR);
                    }
                    
                    // now it's time to try add this guy to the array...
                    // we need to find a suitable kmer...
                    int positioning_kmer_index = 0;
                    bool found_kmer = false;
                    // first find the first non -1 entry, there must be at least CRASS_DEF_NUM_KMERS_4_MODE_HALF
                    while(positioning_kmer_index < numMers4Mode)
                    {
                        if(-1 != (kmer_positions_DR)[positioning_kmer_index])
                            break;
                        positioning_kmer_index++;
                    }
                    // start here, now we look to the differences between this kmer and the next kmer
                    // and make sure that that difference is upheld in the master
                    while(positioning_kmer_index < (numMers4Mode - 1))
                    {
                        if(((kmer_positions_DR)[positioning_kmer_index] - (kmer_positions_DR)[positioning_kmer_index+1]) == ((kmer_positions_DR_master)[positioning_kmer_index] - (kmer_positions_DR_master)[positioning_kmer_index+1]))
                        {
                            found_kmer = true;
                            break;
                        }
                        positioning_kmer_index++;
                    }
                    
                    if(found_kmer)
                    {

                        // note the position of this DR in the array
                        (*DR_offset_map)[*dr_iter] = ((kmer_positions_ARRAY)[positioning_kmer_index] - (kmer_positions_DR)[positioning_kmer_index]);

                        // We need to check that at least CRASS_DEF_PERCENT_IN_ZONE_CUT_OFF percent of bases agree within the "Zone"
                        int this_DR_start_index = 0;
                        int zone_start_index = *dr_zone_start;
                        int comparison_length = (int)tmp_DR.length();
                        
                        // we only need to compare "within" the zone
                        if((*DR_offset_map)[*dr_iter] < *dr_zone_start)
                        {
                            this_DR_start_index = *dr_zone_start - (*DR_offset_map)[*dr_iter];
                        }
                        else if((*DR_offset_map)[*dr_iter] > *dr_zone_start)
                        {
                            zone_start_index = (*DR_offset_map)[*dr_iter];
                        }
                        
                        // work out the comparison length
                        int eff_zone_length = *dr_zone_end - zone_start_index;
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
                            int positional_offset = (kmer_positions_DR_master)[0] - (kmer_positions_DR_master)[positioning_kmer_index] + (kmer_positions_ARRAY)[positioning_kmer_index];
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
										int this_read_start_pos = positional_offset - (*read_iter)->startStopsAt(dr_start_index) - (kmer_positions_DR)[0] ;
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
                                read_iter++;
                            }
                        }
                    }
                }
            }
        }
        dr_iter++;
    }	

    // kill the unfounded ones
    dr_iter = (mDR2GIDMap[GID])->begin();
    while (dr_iter != (mDR2GIDMap[GID])->end()) 
    {
    	if(DR_offset_map->find(*dr_iter) != DR_offset_map->end())
    	{
			if((*DR_offset_map)[*dr_iter] == -1)
			{
				clearReadList(mReads[*dr_iter]);
				mReads[*dr_iter] = NULL;
				dr_iter = (mDR2GIDMap[GID])->erase(dr_iter); 
			}
	    	else
	    	{
	    		dr_iter++;
	    	}
    	}
    	else
    	{
    		dr_iter++;
    	}
    }
    return true;
} 
int WorkHorse::getOffsetAgainstMaster(std::string& masterDR, std::string& slaveDR)
{
    // call smith-waterman
    int a_start_align, a_end_align;
    stringPair pair = smithWaterman(masterDR, slaveDR, &a_start_align, &a_end_align, 0, (int)masterDR.size());
    std::cout<<"a_start_align: "<<a_start_align<< " a_end_align: "<<a_end_align<<std::endl;
    std::cout<<pair.first<<std::endl<<pair.second<<std::endl;
    return a_start_align;
}

std::string WorkHorse::calculateDRConsensus(int GID, std::map<StringToken, int> * DR_offset_map, int * collapsed_pos, std::map<char, int> * collapsed_options, std::map<int, bool> * refined_DR_ends, int * dr_zone_start, int * dr_zone_end, int ** coverage_array, char * consensus_array, float * conservation_array, int * nextFreeGID)
{
	//-----
	// calculate the consensus sequence in the consensus array and the  
	// sequence of the true DR
	// warning, A heavy!
	logInfo("Calculating consensus sequence from aligned reads", 1)
#ifdef DEBUG
		logInfo("DR zone: " << *dr_zone_start << " -> " << *dr_zone_end, 1);
#endif
	
	int array_len = ((int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength > CRASS_DEF_MIN_CONS_ARRAY_LEN) ? (int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength : CRASS_DEF_MIN_CONS_ARRAY_LEN;
	// chars we luv!
    char alphabet[4] = {'A', 'C', 'G', 'T'};
    std::map<char, int> reverse_alphabet;
	for(int i = 0; i < 4; i++)
	{
		reverse_alphabet[alphabet[i]] = i;
	}
	
	// populate the conservation array
    int num_GT_zero = 0;
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
		{
			conservation_array[j] = (float)(max_count)/total_count;
			num_GT_zero++;
		}
		else
		{
			conservation_array[j] = 0;
		}
	}
    
    // check to see that this DR is supported by a bare minimum of reads
    if(num_GT_zero < CRASS_DEF_MIN_READ_DEPTH)
    	logWarn("**WARNING: low confidence DR", 1);
	
	// trim these back a bit (if we trim too much we'll get it back right now anywho)
	*dr_zone_start += CRASS_DEF_DR_ZONE_TRIM_AMOUNT;
	*dr_zone_end -= CRASS_DEF_DR_ZONE_TRIM_AMOUNT;

	// now use this information to find the true direct repeat
	// first work to the left
	while(*dr_zone_start > 0)
	{
		if(conservation_array[(*dr_zone_start) - 1] >= CRASS_DEF_ZONE_EXT_CONS_CUT_OFF)
			(*dr_zone_start)--;
		else
			break;
	}

	// next work to the right
	while(*dr_zone_end < array_len - 1)
	{
		if(conservation_array[(*dr_zone_end) + 1] >= CRASS_DEF_ZONE_EXT_CONS_CUT_OFF)
			(*dr_zone_end)++;
		else
			break;
	}

#ifdef DEBUG
		logInfo("DR zone (post fix): " << *dr_zone_start << " -> " << *dr_zone_end, 1);
#endif
	
	// finally, make the true DR and check for consistency
	std::string true_DR = "";
	
	for(int i = *dr_zone_start; i <= *dr_zone_end; i++)
	{
#ifdef DEBUG
		logInfo("Pos: " << i << " coverage: " << coverage_array[reverse_alphabet[consensus_array[i]]][i] << " conserved(%): " << conservation_array[i] << " consensus: " << consensus_array[i], 1);
#endif		
		(*collapsed_pos)++;
		if(conservation_array[i] >= CRASS_DEF_COLLAPSED_CONS_CUT_OFF)
		{
			(*refined_DR_ends)[i] = true;
			true_DR += consensus_array[i];
		}
		else
		{
			// possible collapsed cluster
			(*refined_DR_ends)[i] = false;
			
	#ifdef DEBUG
			logInfo("-------------", 5); 
			logInfo("Possible collapsed cluster at position: " << *collapsed_pos << " (" << (*dr_zone_start + *collapsed_pos) << " || " << conservation_array[i] << ")", 5);
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
					(*collapsed_options)[alphabet[k]] = (int)(collapsed_options->size() + *nextFreeGID);
					(*nextFreeGID)++;
				}
			}
			
			// make sure we've got more than 1 option
			if(2 > collapsed_options->size())
			{
				collapsed_options->clear();
	#ifdef DEBUG
				logInfo("   ...ignoring (FA)", 5);
	#endif
				true_DR += consensus_array[i];
				(*refined_DR_ends)[i] = true;
			}
			else
			{
				// is this seen at the DR level?
				(*refined_DR_ends)[i] = false;
				std::map<char, int> collapsed_options2;
				DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
				while (dr_iter != (mDR2GIDMap[GID])->end()) 
				{
					std::string tmp_DR = mStringCheck.getString(*dr_iter);
					if(-1 != (*DR_offset_map)[*dr_iter])
					{
						// check if the deciding character is within range of this DR
						// collapsed_pos + dr_zone_start gives the index in the ARRAY of the collapsed char
						// DR_offset_map[*dr_iter] gives the start of the DR in the array
						// We need to check that collapsed_pos + dr_zone_start >= DR_offset_map[*dr_iter] AND
						// that collapsed_pos < dr_zone_start - DR_offset_map[*dr_iter] + tmp_DR.length()
						//if(DR_offset_map[*dr_iter] <= dr_zone_start && dr_zone_start < (DR_offset_map[*dr_iter] + (int)(tmp_DR.length())) && collapsed_pos < (int)(tmp_DR.length()))
						if((*collapsed_pos + *dr_zone_start >= (*DR_offset_map)[*dr_iter]) && (*collapsed_pos + *dr_zone_start - (*DR_offset_map)[*dr_iter] < ((int)tmp_DR.length())))
						{
							// this is easy, we can compare based on this char only
							char decision_char = tmp_DR[*dr_zone_start - (*DR_offset_map)[*dr_iter] + *collapsed_pos];
							collapsed_options2[decision_char] = (*collapsed_options)[decision_char];
						}
					}
					else
					{
						logWarn("No offset for DR: " << tmp_DR, 1);
					}
					dr_iter++;
				}
				
				if(2 > collapsed_options2.size())
				{
					// in the case that the DR is collapsing at the very end of the zone,
					// it may be because the spacers ahve a weird distribution of starting
					// bases. We need to check this out here...
#ifdef DEBUG
					if(*collapsed_pos == 0)
					{
						logInfo("   ...ignoring (RLO SS)", 5);
					}
					else if(*collapsed_pos + *dr_zone_start == *dr_zone_end)
					{
						logInfo("   ...ignoring (RLO EE)", 5);
					}
					else
					{
						logInfo("   ...ignoring (RLO KK)", 5);
					}
#endif
					true_DR += consensus_array[i];
					(*refined_DR_ends)[i] = true;
					collapsed_options->clear();
				}
				else
				{
					// If it aint in collapsed_options2 it aint super fantastic!
					collapsed_options->clear();
					*collapsed_options = collapsed_options2;
					
					// make the collapsed pos array specific and exit this loop
					*collapsed_pos += *dr_zone_start;
					i = *dr_zone_end + 1;
				}
			}
		}
	}
	logInfo("Consensus DR: " << true_DR, 1);
	return true_DR;
}

bool WorkHorse::parseGroupedDRs(int numMers4Mode, int GID, std::vector<std::string> * nTopKmers, int * nextFreeGID) 
{
	
    //-----
    // Cluster refinement and possible splitting for a Group ID
    //
    logInfo("Parsing group: " << GID, 4);
    		
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Find a Master DR for this group of DRs
    StringToken master_DR_token = -1;
    std::string master_DR_sequence = "**unset**";
    if(!findMasterDR(GID, nTopKmers, &master_DR_token, &master_DR_sequence)) { return false; }
    
    // now we have the n most abundant kmers and one DR which contains them all
    // time to rock and rrrroll!

    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialise variables we'll need
    // chars we luv!
    char alphabet[4] = {'A', 'C', 'G', 'T'};
    int array_len = ((int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength > CRASS_DEF_MIN_CONS_ARRAY_LEN) ? (int)CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mAveReadLength : CRASS_DEF_MIN_CONS_ARRAY_LEN;

    // first we need a 4 * array_len
    int ** coverage_array = new int*[4];
    
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
    
    // The offset of the start position of each potential DR 
    // when compared to the "true DR"
    // we use this structure when we detect overcollapsing
    std::map<StringToken, int> DR_offset_map;
    
    // we will compare the orientations and positions of all other guys to the master so we need a permanent store
    bool kmer_rcs_DR_master[CRASS_DEF_NUM_KMERS_4_MODE];
    int kmer_positions_DR_master[CRASS_DEF_NUM_KMERS_4_MODE];
    
    // look for the start and end of the DR zone
    int dr_zone_start = -1;
    int dr_zone_end = -1;

    // just the positions in the DR fangs...
    kmer_positions_ARRAY[0] = (int)(array_len*CRASS_DEF_CONS_ARRAY_START);
    isKmerPresent(kmer_rcs_DR, kmer_positions_DR, (*nTopKmers)[0], &master_DR_sequence);
    
    for(int i = 1; i < numMers4Mode; i++)
    {
    	//std::cout << "1" << (kmer_rcs_DR + i) << std::endl;
    	//std::cout << "2" << (kmer_positions_DR + i) << std::endl;
    	//std::cout << "3" << (*nTopKmers)[i] << std::endl;
    	//std::cout << "4" << &master_DR_sequence << std::endl;
        isKmerPresent((kmer_rcs_DR + i), (kmer_positions_DR + i), (*nTopKmers)[i], &master_DR_sequence);
        kmer_positions_ARRAY[i] = kmer_positions_DR[i] - kmer_positions_DR[0] + kmer_positions_ARRAY[0];
    }
    
    // store the first results away as the master results
    for(int i = 0; i < numMers4Mode; i++)
    {
        kmer_rcs_DR_master[i] = kmer_rcs_DR[i];
        kmer_positions_DR_master[i] = kmer_positions_DR[i];
    }

    // note the position of the master DR in the array
    DR_offset_map[master_DR_token] = (kmer_positions_ARRAY[0] - kmer_positions_DR[0]);

    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Set up the master DR's array and insert this guy into the main array
    populateCoverageArray(numMers4Mode, GID, master_DR_sequence, master_DR_token, &DR_offset_map, &dr_zone_start, &dr_zone_end, nTopKmers, coverage_array, kmer_positions_DR, kmer_rcs_DR, kmer_positions_DR_master, kmer_rcs_DR_master, kmer_positions_ARRAY);
    
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // calculate consensus and diversity
	// use these variables to identify and store possible
	// collapsed clusters
	int collapsed_pos = -1;
	std::map<char, int> collapsed_options;            // holds the chars we need to split on
	std::map<int, bool> refined_DR_ends;              // so we can update DR ends based on consensus 
    std::string true_DR = calculateDRConsensus(GID, &DR_offset_map, &collapsed_pos, &collapsed_options, &refined_DR_ends, &dr_zone_start, &dr_zone_end, coverage_array, consensus_array, conservation_array, nextFreeGID);

    // check to make sure that the DR is not just some random long RE
    if((unsigned int)(true_DR.length()) > mOpts->highDRsize)
    {
        // probably a dud. throw it out
        // free the memory and clean up
        if(NULL != mDR2GIDMap[GID])
        {
        	delete mDR2GIDMap[GID];
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
        if(NULL != mDR2GIDMap[GID])
        {
        	delete mDR2GIDMap[GID];
        	mDR2GIDMap[GID] = NULL;
        }
        
        logInfo("Killed: {" << true_DR << "} cause' the consensus was too short... (" << true_DR.length() << " ," << collapsed_options.size() << ")", 1);

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
    
    
    // QC the DR again for low complexity
    if (isRepeatLowComplexity(true_DR) && collapsed_options.size() == 0) 
    {
        // probably a dud. throw it out
        // free the memory and clean up
        if(NULL != mDR2GIDMap[GID])
        {
        	delete mDR2GIDMap[GID];
        	mDR2GIDMap[GID] = NULL;
        }
        
        logInfo("Killed: {" << true_DR << "} cause' the consensus was low complexity...", 1);
        
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
    // Update the refined starts and ends of the DR 
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
        
        // print out the consensus array
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
        logInfo("Attempting to split the collapsed DR", 5);
        std::map<char, int> coll_char_to_GID_map;
        std::map<char, int>::iterator co_iter = collapsed_options.begin();
        while(co_iter != collapsed_options.end())
        {
            int group = (*nextFreeGID)++;
            mDR2GIDMap[group] = new DR_Cluster;
            coll_char_to_GID_map[co_iter->first] = group;
            logInfo("Mapping \""<< co_iter->first << " : "  << co_iter->second << "\" to group: " << group, 1);
            co_iter++;
        }
        
        DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
        while (dr_iter != (mDR2GIDMap[GID])->end()) 
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
        if(NULL != mDR2GIDMap[GID])
        {
        	delete mDR2GIDMap[GID];
        	mDR2GIDMap[GID] = NULL;
        }
        
        logInfo("Calling the parser recursively", 4);
        
        // call this baby recursively with the new clusters
        std::map<char, int>::iterator cc_iter = coll_char_to_GID_map.begin();
        while(cc_iter != coll_char_to_GID_map.end())
        {
            parseGroupedDRs(numMers4Mode, cc_iter->second, nTopKmers, nextFreeGID);
            cc_iter++;
        }
    }
    else
    {
        //++++++++++++++++++++++++++++++++++++++++++++++++
        // repair all the startstops for each read in this group
        //
        // This function is recursive, so we'll only get here when we have found exactly one DR
        
        // make sure that the true DR is in its laurenized form
        std::string laurenized_true_dr = laurenize(true_DR);
        bool rev_comp = (laurenized_true_dr != true_DR) ? true : false;

        logInfo("Found DR: " << laurenized_true_dr, 2);
        
        mTrueDRs[GID] = laurenized_true_dr;
        DR_ClusterIterator drc_iter = (mDR2GIDMap[GID])->begin();
        while(drc_iter != (mDR2GIDMap[GID])->end())
        {
        	if(DR_offset_map.find(*drc_iter) == DR_offset_map.end())
        	{
        		logError("1: Repeat "<< *drc_iter<<" in Group "<<GID <<" has no offset in DR_offset_map");
        	}
        	else
        	{
				if (DR_offset_map[*drc_iter] == -1 ) 
				{
					// This means that we couldn't add this DR or any of it's reads into the consensus array
					logError("2: Repeat "<< *drc_iter<<" in Group "<<GID <<" has no offset in DR_offset_map");
				} 
				else 
				{
					// go through each read
					ReadListIterator read_iter = mReads[*drc_iter]->begin();
					while (read_iter != mReads[*drc_iter]->end()) 
					{
						(*read_iter)->updateStartStops((DR_offset_map[*drc_iter] - dr_zone_start), &true_DR, mOpts);
	
						// reverse complement sequence if the true DR is not in its laurenized form
						if (rev_comp) 
						{
							(*read_iter)->reverseComplementSeq();
						}
						read_iter++;
					}
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

bool WorkHorse::isKmerPresent(bool * didRevComp, int * startPosition, const std::string kmer, const std::string *  DR)
{
    //-----
    // Work out if a Kmer is present in a string and store positions etc...
    //
    std::string tmp_kmer = reverseComplement(kmer);
    size_t pos = DR->find(kmer);
    if(pos == string::npos)
    {
        // try the reverse complement
        // rev compt the kmer, it's shorter!
        pos = DR->find(tmp_kmer);
        if(pos != string::npos)
        {
            // found the kmer in the reverse direction!
        	// make sure I found it once only
            size_t pos2 = DR->find(tmp_kmer, pos+1);
            if(pos2 != string::npos)
            {
            	// found it twice
                *startPosition = -1;
                return false;
            } // else OK
            
            *didRevComp = true;
            *startPosition = (int)pos;    
            return true;
        }
    }
    else
    {
        // found the kmer in the forward direction!
        // search for more in the forward direction
        size_t pos2 = DR->find(kmer, pos+1);
        if(pos2 != string::npos)
        {
        	// found it twice    	
            *startPosition = -1;
            return false;
        }
        else
        {
            // none? -> search in the reverse direction from start
            pos2 = DR->find(tmp_kmer);
            if(pos2 != string::npos)
            {
            	// found in both dirs!
                *startPosition = -1;
                return false;
            } // else OK
        }
        
        *didRevComp = false;
        *startPosition = (int)pos;
        return true;
    }
    *startPosition = -1;
    return false;
}

int WorkHorse::getNMostAbundantKmers(std::vector<std::string>& mostAbundantKmers, int num2Get, std::map<std::string, int> * kmer_CountMap)
{
	//-----
	// True to it's name get MOST abundant kmers.
	// 
	return getNMostAbundantKmers(1000000, mostAbundantKmers, num2Get, kmer_CountMap);
}

int WorkHorse::getNMostAbundantKmers(int maxAmount, std::vector<std::string>& mostAbundantKmers, int num2Get, std::map<std::string, int> * kmer_CountMap)
{
	//-----
	// get the most abundant kmers under a certain amount.
	//
    std::string top_kmer;    
    std::map<std::string, bool> top_kmer_map;
    
    if ((int)(kmer_CountMap->size()) < num2Get) 
    {
        return 0;
    } 
    else 
    {
        for (int i = 1; i <= num2Get; i++) 
        {
            std::map<std::string, int>::iterator map_iter = kmer_CountMap->begin();
            int max_count = 0;
            
            while (map_iter != kmer_CountMap->end()) 
            {
            	//std::cout << map_iter->first << " : " << map_iter->second << " : " << max_count << " : " << maxAmount << std::endl;
                if((map_iter->second > max_count) && (map_iter->second <= maxAmount) && (top_kmer_map.find(map_iter->first) == top_kmer_map.end()))
                {
                    max_count = map_iter->second;
                    top_kmer = map_iter->first;
                    //std::cout << "NT: " << top_kmer << std::endl;
                }
                map_iter++;
            }
            //std::cout << "ADDING: " << top_kmer << std::endl;            
            top_kmer_map[top_kmer] = true;
        }
        int num_mers_found = 0;
        std::map<std::string, bool>::iterator tkm_iter = top_kmer_map.begin();
        while(tkm_iter != top_kmer_map.end())
        {
            //std::cout<<tkm_iter->first<<std::endl;
        	num_mers_found++;
            mostAbundantKmers.push_back(tkm_iter->first);
            tkm_iter++;
        }
        return num_mers_found;
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
    int min_clust_membership_count = mOpts->kmer_clust_size;
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
        
        // see if this guy has been counted LOCALLY
        // use this list when we go to select N most abundant kmers
        if(local_kmer_CountMap.find(tmp_str) == local_kmer_CountMap.end())
        {
        	// first time this kmer has been seen in this read
            local_kmer_CountMap[tmp_str] = 1;
        }
        else
        {
            local_kmer_CountMap[tmp_str]++;
        }
        
        // see if we've seen this kmer before GLOBALLY
        std::map<std::string, int>::iterator k2g_iter = k2GIDMap->find(tmp_str);
        if(k2g_iter == k2GIDMap->end())
        {
            // first time we seen this one GLOBALLY
            homeless_kmers.push_back(tmp_str);
        }
        else
        {
        	// we've seen this guy before.
            // only do this if our guy doesn't belong to a group yet
            if(0 == group)
            {
                // this kmer belongs to a group -> increment the local group count
                std::map<int, int>::iterator this_group_iter = group_count.find(k2g_iter->second);
                if(this_group_iter == group_count.end())
                {
                    group_count[k2g_iter->second] = 1;
                }
                else
                {
                    group_count[k2g_iter->second]++;
                    // have we seen this guy enought times?
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
    	// we couldn't put our guy into a group
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
                logInfo("Assigning flankers for NodeManager "<<drg_iter->first, 3);
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

//int WorkHorse::dumpSpacers(void)
//{
//	//-----
//	// Wrapper for graph cleaning
//	//
//	// create a spacer dictionary
//	logInfo("Dumping spacers", 1);
//	DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
//	while(drg_iter != mDR2GIDMap.end())
//	{
//		if(NULL != drg_iter->second)
//		{            
//			if (NULL != mDRs[mTrueDRs[drg_iter->first]])
//            {
//                (mDRs[mTrueDRs[drg_iter->first]])->dumpSpacerDict(mOpts->output_fastq + "Group_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first] + ".spacers", false);
//		    }
//        }
//		drg_iter++;
//	}
//	return 0;
//}

//void WorkHorse::writeLookupToFile(string &outFileName, lookupTable &outLookup)
//{
//    std::ofstream outLookupFile;
//    outLookupFile.open(outFileName.c_str());
//    
//    lookupTable::iterator ter = outLookup.begin();
//    while (ter != outLookup.end()) 
//    {
//        outLookupFile<<ter->first<<"\t"<<ter->second<<endl;
//        
//        ter++;
//    }
//    outLookupFile.close();
//}

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
                
                // check to see if there is anything to print - if so then make the key
                if ( mDRs[mTrueDRs[drg_iter->first]]->printSpacerGraph(graph_file_name, mTrueDRs[drg_iter->first], mOpts->longDescription, mOpts->showSingles))
                {
                    mDRs[mTrueDRs[drg_iter->first]]->printSpacerKey(key_file, 10, namePrefix + to_string(drg_iter->first));
                }
                else 
                {
                    // should delete this guy since there are no spacers
                    // this way the group will not be in the xml either
                    delete mDRs[mTrueDRs[drg_iter->first]];
                    mDRs[mTrueDRs[drg_iter->first]] = NULL;
                }
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
	
    crispr::XML * xml_doc = new crispr::XML();
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

bool WorkHorse::addDataToDOM(crispr::XML * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    try 
    {
        xercesc::DOMElement * data_elem = xmlDoc->addData(groupElement);
        if ((mDRs[mTrueDRs[groupNumber]])->haveAnyFlankers()) {
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

bool WorkHorse::addMetadataToDOM(crispr::XML * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    try{
        
        std::stringstream notes;
        notes << PACKAGE_NAME <<" ("<<PACKAGE_VERSION<<") run on "<<mTimeStamp<<" with command: ";
        notes <<mCommandLine;
        xercesc::DOMElement * metadata_elem = xmlDoc->addMetaData(notes.str(), groupElement);
        
        std::string file_name;
        char * buf = NULL;
        std::string absolute_dir = getcwd(buf, 4096);
        absolute_dir += "/";
        delete buf;
        // add in files if they exist
        if (!mOpts->logToScreen) 
        {
            // we whould have a log file
            file_name = mOpts->output_fastq + PACKAGE_NAME + "." + mTimeStamp + ".log";
            if (checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("log", absolute_dir + file_name, metadata_elem);
            }
            else
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
        }
        
        
    #ifdef DEBUG
        // check for debuging .gv files
        if (!mOpts->noDebugGraph) 
        {
            file_name = mOpts->output_fastq + "Group_"; 
            std::string file_sufix = to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + "_debug.gv";
            if (checkFileOrError((file_name + file_sufix).c_str())) 
            {
                xmlDoc->addFileToMetadata("data", absolute_dir + file_name + file_sufix, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, (absolute_dir + file_name + file_sufix).c_str() );
            }
            
            // and now for the cleaned .gv
            file_name = mOpts->output_fastq + "Clean_";
            if (checkFileOrError((file_name + file_sufix).c_str())) 
            {
                xmlDoc->addFileToMetadata("data", absolute_dir + file_name + file_sufix, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, (absolute_dir + file_name + file_sufix).c_str() );                
            }
        }

        
    #endif // DEBUG
        
    #ifdef RENDERING
        // check for image files
    #ifdef DEBUG
        if (!mOpts->noDebugGraph) 
        {
            file_name = mOpts->output_fastq + "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
            if (checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("image", absolute_dir + file_name, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
            
            file_name = mOpts->output_fastq + "Clean_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
            
            if (checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("image", absolute_dir + file_name, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
        }

    #endif // DEBUG
        if (!mOpts->noRendering) 
        {
            file_name = mOpts->output_fastq + "Spacers_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
            if (checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("image", absolute_dir + file_name, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
        } 


    #endif // RENDERING
        
        // check the sequence file
        file_name = mOpts->output_fastq +  "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".fa";
        if (checkFileOrError(file_name.c_str())) 
        {
            xmlDoc->addFileToMetadata("sequence", absolute_dir + file_name, metadata_elem);
        } 
        else 
        {
            throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
        }
    } catch(crispr::no_file_exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    } catch(std::exception& e) {
        std::cerr<<e.what()<<std::endl;
        return 1;
    } catch(xercesc::DOMException& e) {
        
    }
    return 0;
    
}



