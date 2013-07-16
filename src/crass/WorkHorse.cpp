// File: WorkHorse.cpp
// Original Author: Michael Imelfort 2011
// Hacked to death by Connor Skennerton 2012
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Implementation of WorkHorse functions
// 
// --------------------------------------------------------------------
//  Copyright  2011, 2012 Michael Imelfort and Connor Skennerton
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
#include <ctime>
#include <libcrispr/StlExt.h>
#include <libcrispr/Exception.h>

// local includes
#include "WorkHorse.h"
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crassDefines.h"
#include "NodeManager.h"
#include "ReadHolder.h"
#include "SeqUtils.h"
#include "SmithWaterman.h"
#include "StringCheck.h"
#include "config.h"
#include "ksw.h"

bool sortLengthDecending( const std::string& a, const std::string& b)
{
    return a.length() > b.length();
}

bool sortLengthAssending( const std::string& a, const std::string& b)
{
    return a.length() < b.length();
}

// a should be shorter than b if sorted correctly
bool includeSubstring(const std::string& a, const std::string& b)
{
    if (std::string::npos != b.find(a)) {
        return true;
    } else if(std::string::npos != b.find(reverseComplement(a))) {
        return true;
    }
    return false;
}

bool isNotEmpty(const std::string& a)
{
    return !a.empty();
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
int WorkHorse::doWork(Vecstr seqFiles)
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
#ifdef SEARCH_SINGLETON
    std::ofstream debug_out;
    std::stringstream debug_out_file_name;
    debug_out_file_name <<mOpts->output_fastq<< "crass.debug."<<mTimeStamp<<".report";
    debug_out.open((debug_out_file_name.str()).c_str());
    if(debug_out.good()) {
        SearchCheckerList::iterator debug_iter;
        for (debug_iter = debugger->begin(); debug_iter != debugger->end(); debug_iter++) {
            debug_out<<debug_iter->first<<"\t"<<debug_iter->second.gid()<<"\t"<<debug_iter->second.truedr()<<"\t";
            std::vector<StringToken>::iterator node_iter = debug_iter->second.begin();
            if(node_iter != debug_iter->second.end()) {
                debug_out <<*node_iter;
                for ( node_iter++ ; node_iter != debug_iter->second.end(); node_iter++) {
                    debug_out<<":"<<*node_iter;
                }
            }
            debug_out<<"\t";
            Vecstr::iterator sp_iter = debug_iter->second.beginSp();
            if(sp_iter != debug_iter->second.endSp()) {
                debug_out<<*sp_iter;
                for ( sp_iter++ ; sp_iter != debug_iter->second.endSp(); sp_iter++) {
                    debug_out<<":"<<*sp_iter;
                }
            }
            debug_out<<std::endl;
        }
    } else {
        std::cerr<<"error printing debugging report"<<std::endl;
        return 200;
    }
#endif
	
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
    // call flanking regions
    if (generateFlankers()) {
        logError("FATAL ERROR: generateFlankers failed");
        return 70;
    }
    
    //remove NodeManagers with low numbers of spacers
    // and where the standard deviation of the spacer length 
    // is too high
    if (removeLowConfidenceNodeManagers())
    {
        logError("FATAL ERROR: removeLowSpacerNodeManagers failed");
        return 7;
    }
	
    // print the reads to a file if requested
//	if(dumpReads(false))
//	{
//        logError("FATAL ERROR: dumpReads failed");
//        return 9;
//	}
	
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
//	if(renderSpacerGraphs())
//	{
//        logError("FATAL ERROR: renderSpacerGraphs failed");
//        return 11;
//	}
	
	outputResults();
	
    logInfo("all done!", 1);
	return 0;
}

int WorkHorse::parseSeqFiles(Vecstr seqFiles)
{
	//-----
	// Load data from files and search for DRs
	//
    Vecstr::iterator seq_iter = seqFiles.begin();
    
    // direct repeat sequence and unique ID
    lookupTable patterns_lookup;
    
    // the sequence of whole spacers and their unique ID
    lookupTable reads_found;

    time_t start_time;
    time(&start_time);
    while(seq_iter != seqFiles.end())
    {
        logInfo("Parsing file: " << *seq_iter, 1);
        try {
            int max_len = decideWhichSearch(seq_iter->c_str(), 
                                            *mOpts, 
                                            &mReads, 
                                            &mStringCheck, 
                                            patterns_lookup, 
                                            reads_found,
                                            start_time);
            
            mMaxReadLength = (max_len > mMaxReadLength) ? max_len : mMaxReadLength;
            logInfo("Finished file: " << *seq_iter, 1);

        } catch (crispr::exception& e) {
            std::cerr<<e.what()<<std::endl;
            return 1;
        }
        
        // Check to see if we found anything, should return if we haven't
        if (patterns_lookup.empty()) 
        {
            logInfo("No direct repeat sequences were identified for file: "<<seq_iter->c_str(), 1);
        }
        logInfo("Finished file: " << *seq_iter, 1);
        
        seq_iter++;
    }
    // add in a new line so the looger won't overlap itself
    std::cout<<std::endl;

    GroupKmerMap group_kmer_counts_map;
    int next_free_GID = 1;
    Vecstr * non_redundant_set = createNonRedundantSet(group_kmer_counts_map, next_free_GID);
    logInfo("Number of reads found so far: "<<this->numOfReads(), 2);

    if (non_redundant_set->size() > 0) 
    {
        std::cout<<"["<<PACKAGE_NAME<<"_clusterCore]: " << non_redundant_set->size() << " non-redundant patterns."<<std::endl;
        seq_iter = seqFiles.begin();
        logInfo("Begining Second iteration through files to recruit singletons", 2);


        time(&start_time);
        while (seq_iter != seqFiles.end()) {
            
            logInfo("Parsing file: " << *seq_iter, 1);
            
            try {
                findSingletons(seq_iter->c_str(), *mOpts, non_redundant_set, reads_found, &mReads, &mStringCheck, start_time);
            } catch (crispr::exception& e) {
                std::cerr<<e.what()<<std::endl;
                delete non_redundant_set;
                return 1;
            }
            seq_iter++;
        }
    }
    // add in a new line so the ouptut won't overlap itself
    std::cout<<std::endl;
    delete non_redundant_set;
    std::cout<<"["<<PACKAGE_NAME<<"_patternFinder]: "<<"Found "<<numOfReads()<<" reads"<<std::endl;
    logInfo("Searching complete. " << mReads.size()<<" direct repeat variants have been found", 1);
    logInfo("Number of reads found so far: "<<this->numOfReads(), 2);

    if(mOpts->removeHomopolymers) {
        // change back the sizes of the direct repeats to 
        // counter the changes made by mOpts.removeHomopolymers
        // this way the final DRs and spacer should fall 
        // inside the correct lengths
        mOpts->lowDRsize /= mOpts->averageDrScalling;
        mOpts->highDRsize /= mOpts->averageDrScalling;
        mOpts->lowSpacerSize /= mOpts->averageSpacerScalling;
        mOpts->highSpacerSize /= mOpts->averageSpacerScalling;
    }
    try {
        if (findConsensusDRs(group_kmer_counts_map, next_free_GID))
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
    //std::cout<<'['<<PACKAGE_NAME<<"_graphBuilder]: "<<mTrueDRs.size()<<" putative CRISPRs found!"<<std::endl;
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
                    if(*read_iter == NULL) {
                        logError("Read is set to null");
                    }
                    //MI std::cout<<'.'<<std::flush;
#ifdef SEARCH_SINGLETON
                    SearchCheckerList::iterator debug_iter = debugger->find((*read_iter)->getHeader());
                    if (debug_iter != debugger->end()) {
                        //found one of our interesting reads
                        // add in the true DR
                        debug_iter->second.truedr(mTrueDRs[drg_iter->first]);
                        debug_iter->second.gid(drg_iter->first);
                    }
#endif
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

int WorkHorse::removeLowConfidenceNodeManagers(void)
{
    logInfo("Removing CRISPRs with low numbers of spacers", 1);
	int counter = 0;
    DR_Cluster_MapIterator drg_iter = mDR2GIDMap.begin();
	while(drg_iter != mDR2GIDMap.end())
	{
		if(NULL != drg_iter->second)
		{            
            if (NULL != mDRs[mTrueDRs[drg_iter->first]])
            {
                NodeManager * current_manager = mDRs[mTrueDRs[drg_iter->first]];
                if( current_manager->getSpacerCountAndStats(false) < mOpts->covCutoff) 
                {
                    logInfo("Deleting NodeManager "<<drg_iter->first<<" as it contained less than "<<mOpts->covCutoff<<" attached spacers",5);
                    delete mDRs[mTrueDRs[drg_iter->first]];
                     mDRs[mTrueDRs[drg_iter->first]] = NULL;
                } else if (current_manager->stdevSpacerLength() > CRASS_DEF_STDEV_SPACER_LENGTH) {
                    logInfo("Deleting NodeManager "<<drg_iter->first<<" as the stdev ("<<current_manager->stdevSpacerLength()<<") of the spacer lengths was greater than "<<CRASS_DEF_STDEV_SPACER_LENGTH, 4);
                    delete mDRs[mTrueDRs[drg_iter->first]];
                     mDRs[mTrueDRs[drg_iter->first]] = NULL;
                }
                counter++;
            }
		}
		drg_iter++;
	}
    std::cout<<'['<<PACKAGE_NAME<<"_graphBuilder]: "<<counter<<" putative CRISPRs have passed all checks"<<std::endl;
	return 0;
}

//**************************************
// Functions used to cluster DRs into groups and identify the "true" DR
//**************************************
int WorkHorse::findConsensusDRs(GroupKmerMap& groupKmerCountsMap, int& nextFreeGID)
{
    //-----
    // Cluster potential DRs and work out their true sequences
    // make the node managers while we're at it!
    //

    logInfo("Reducing list of potential DRs (2): Cluster refinement and true DR finding", 1);
    
    // go through all the counts for each group
    GroupKmerMap::iterator group_count_iter; 
    for(group_count_iter =  groupKmerCountsMap.begin(); 
        group_count_iter != groupKmerCountsMap.end(); 
        group_count_iter++)
    {
        if(NULL == mDR2GIDMap[group_count_iter->first])
        {
            continue;
        }

        parseGroupedDRs(group_count_iter->first, &nextFreeGID);
        
        // delete the kmer count lists cause we're finsihed with them now
        if(NULL != group_count_iter->second)
        {
            delete group_count_iter->second;
            group_count_iter->second = NULL;
        }
    }
    
    return 0;
}
void WorkHorse::removeRedundantRepeats(Vecstr& repeatVector)
{
    // given a vector of repeat sequences, will order the vector based on repeat
    // length and then remove longer repeats if there is a shorter one that is
    // a perfect substring
    std::sort(repeatVector.begin(), repeatVector.end(), sortLengthAssending);
    Vecstr::iterator iter;

    // go though all of the patterns and determine which are substrings
    // clear the string if it is
    for (iter = repeatVector.begin(); iter != repeatVector.end(); iter++) {
        if (iter->empty()) {
            continue;
        }
        Vecstr::iterator iter2;
        for (iter2 = iter+1; iter2 != repeatVector.end(); iter2++) {
            if (iter2->empty()) {
                continue;
            }
            //pass both itererators into the comparison function
            if (includeSubstring(*iter, *iter2)) {
                iter2->clear();
            }
        }

    }

    // ok so now partition the vector so that all the empties are at one end
    // will return an iterator postion to the first position where the string
    // is empty
    Vecstr::iterator empty_iter = std::partition(repeatVector.begin(), repeatVector.end(), isNotEmpty);
    // remove all the empties from the list
    repeatVector.erase(empty_iter, repeatVector.end());
}


Vecstr * WorkHorse::createNonRedundantSet(GroupKmerMap& groupKmerCountsMap, int& nextFreeGID)
{
    // cluster the direct repeats then remove the redundant ones
    // creates a vector in dynamic memory, so don't forget to delete 
    //-----
    // Cluster potential DRs and work out their true sequences
    // make the node managers while we're at it!
    //
    std::map<std::string, int> k2GID_map;
    logInfo("Reducing list of potential DRs (1): Initial clustering", 1);
    logInfo("Reticulating splines...", 1);    
    // go through all of the read holder objects
    ReadMapIterator read_map_iter = mReads.begin();
    while (read_map_iter != mReads.end()) 
    {
        clusterDRReads(read_map_iter->first, &nextFreeGID, &k2GID_map, &groupKmerCountsMap);
        ++read_map_iter;
    }
    std::cout<<'['<<PACKAGE_NAME<<"_clusterCore]: "<<mReads.size()<<" variants mapped to "<<mDR2GIDMap.size()<<" clusters"<<std::endl;
    std::cout<<'['<<PACKAGE_NAME<<"_clusterCore]: creating non-redundant set"<<std::endl;

    Vecstr * non_redundant_repeats = new Vecstr();

    DR_Cluster_MapIterator dcg_iter = mDR2GIDMap.begin();
    while(dcg_iter != mDR2GIDMap.end())
    {
        DR_ClusterIterator dc_iter = (dcg_iter->second)->begin();
        if (dcg_iter->second != NULL) 
        {
            logInfo("-------------", 4);
            logInfo("Group: " << dcg_iter->first, 4);
            
            Vecstr clustered_repeats;
            while(dc_iter != (dcg_iter->second)->end())
            {
                std::string tmp = mStringCheck.getString(*dc_iter);
                clustered_repeats.push_back(tmp);
                logInfo(tmp, 4);
                dc_iter++;
            }
            logInfo("-------------", 4);
            
            removeRedundantRepeats(clustered_repeats);
            Vecstr::iterator cr_iter;
            Vecstr tmp_vec;
            for (cr_iter = clustered_repeats.begin(); cr_iter != clustered_repeats.end(); ++cr_iter) {
                tmp_vec.push_back(reverseComplement(*cr_iter));
            }
            non_redundant_repeats->insert(non_redundant_repeats->end(), clustered_repeats.begin(), clustered_repeats.end());
            non_redundant_repeats->insert(non_redundant_repeats->end(), tmp_vec.begin(), tmp_vec.end());

        }
        dcg_iter++;
    }
    logInfo("non-redundant patterns:", 4);
    Vecstr::iterator nr_iter;
    for (nr_iter = non_redundant_repeats->begin(); nr_iter != non_redundant_repeats->end(); ++nr_iter) {
        logInfo(*nr_iter, 4);
    }
    logInfo("-------------", 4);
    return non_redundant_repeats;
}

bool WorkHorse::findMasterDR(int GID, StringToken * masterDRToken, std::string * masterDRSequence)
{
	//-----
	// Identify a master DR
	// 
	// Updates the values in masterDRToken and  masterDRSequence and returns true
	// otherwise deletes the memory pointed at by clustered_DRs and returns false
	//
	//
    
    // this new version of the function needs only to find the longest DR in the list
    // so that all of the other DRs can be aligned against it.
    
    logInfo("Identifying a master DR", 1);

    DR_Cluster * current_dr_cluster = mDR2GIDMap[GID];
    size_t current_longest_size = 0;
    DR_ClusterIterator dr_iter;// = dr_cluster->begin();
    
    for (dr_iter = current_dr_cluster->begin(); dr_iter != current_dr_cluster->end(); ++dr_iter) {

        std::string tmp_dr_seq = mStringCheck.getString(*dr_iter);
        if (tmp_dr_seq.size() > current_longest_size) {
            *masterDRToken = *dr_iter;
            *masterDRSequence = tmp_dr_seq;
            current_longest_size = tmp_dr_seq.size();
        }   
    }
    if(*masterDRToken == -1)
    {
        logError("Could not identify a master DR");
    }
    logInfo("Identified: " << *masterDRSequence << " (" << *masterDRToken << ") as a master potential DR", 4);


    
    return true;
}

bool WorkHorse::populateCoverageArray(int GID, Aligner& drAligner)
{
	//-----
	// Use the data structures initialised in parseGroupedDRs
	// Load all the reads into the consensus array
	//
	logInfo("Populating consensus array", 1);

    //++++++++++++++++++++++++++++++++++++++++++++++++
    // now go thru all the other DRs in this group and add them into
    // the consensus array
    DR_ClusterIterator dr_iter;// = (mDR2GIDMap[GID])->begin();
    for (dr_iter = (mDR2GIDMap[GID])->begin(); dr_iter != (mDR2GIDMap[GID])->end(); dr_iter++) 
    {
        // we've already done the master DR
        if(drAligner.getMasterDrToken() == *dr_iter)
        {
            continue;
        }

        // get the string for this mofo
        //std::string tmp_DR = mStringCheck.getString(*dr_iter);
        drAligner.alignSlave(*dr_iter);
    }
    // kill the unfounded ones
    dr_iter = (mDR2GIDMap[GID])->begin();
    while (dr_iter != (mDR2GIDMap[GID])->end()) 
    {
    	if(drAligner.offsetFind(*dr_iter) != drAligner.offsetEnd())
    	{
			if(drAligner.offset(*dr_iter) == -1)
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


std::string WorkHorse::calculateDRConsensus(int GID, 
                                            Aligner& drAligner, 
                                            int& nextFreeGID,
                                            int& collapsedPos,
                                            std::map<char, int>& collapsedOptions,
                                            std::map<int, bool>& refinedDREnds
                                            )
{

    drAligner.generateConsensus();
	
	// make the true DR and check for consistency
	std::string true_dr = "";
    for (int i = drAligner.getDRZoneStart(); i <= drAligner.getDRZoneEnd(); i++)
	{
#ifdef DEBUG
        logInfo("Pos: " << i << " coverage: " << drAligner.coverageAt(i, drAligner.consensusAt(i)) << " conserved(%): " << drAligner.conservationAt(i) << " consensus: " << drAligner.consensusAt(i), 1);
#endif		
		collapsedPos++;
		if(drAligner.conservationAt(i) >= CRASS_DEF_COLLAPSED_CONS_CUT_OFF)
		{
			refinedDREnds[i] = true;
			true_dr += drAligner.consensusAt(i);
		}
		else
		{
			// possible collapsed cluster
			refinedDREnds[i] = false;
			
	#ifdef DEBUG
			logInfo("-------------", 5); 
			logInfo("Possible collapsed cluster at position: " << collapsedPos << " (" << (drAligner.getDRZoneStart() + collapsedPos) << " || " << drAligner.conservationAt(i) << ")", 5);
			logInfo("Base:  Count:  Cov:",5);
	#endif
			float total_count = drAligner.depthAt(i);
			char alphabet[4] = {'A', 'C', 'G', 'T'};
			for(int k = 0; k < 4; k++)
			{
                float nt_proportion = static_cast<float>(drAligner.coverageAt(i, alphabet[k])/total_count);

				// check to make sure that each base is represented enough times
				if(nt_proportion >= CRASS_DEF_COLLAPSED_THRESHOLD)
				{
#ifdef DEBUG
                logInfo("  " << alphabet[k] << "     " << drAligner.coverageAt(i, alphabet[k]) << "      " << nt_proportion <<"   above threshold, added to colapsed options", 5);
#endif
                    // there's enough bases here to warrant further investigation
					collapsedOptions[alphabet[k]] = (int)(collapsedOptions.size() + nextFreeGID);
					nextFreeGID++;
				}
#ifdef DEBUG
                else
                {
                    logInfo("  " << alphabet[k] << "     " << drAligner.coverageAt(i, alphabet[k]) << "      " << nt_proportion <<"   below threshold", 5);
                }
#endif
			}
			
			// make sure we've got more than 1 option
			if(2 > collapsedOptions.size())
			{
				collapsedOptions.clear();
#ifdef DEBUG
				logInfo("   ...ignoring (FA)", 5);
#endif
				true_dr += drAligner.consensusAt(i);
				refinedDREnds[i] = true;
			}
			else
			{
				// is this seen at the DR level?
				refinedDREnds[i] = false;
				std::map<char, int> collapsed_options2;
				DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
				while (dr_iter != (mDR2GIDMap[GID])->end()) 
				{
					std::string tmp_DR = mStringCheck.getString(*dr_iter);
					if(-1 != drAligner.offset(*dr_iter))
					{
						// check if the deciding character is within range of this DR
						// collapsed_pos + dr_zone_start gives the index in the ARRAY of the collapsed char
						// DR_offset_map[*dr_iter] gives the start of the DR in the array
						// We need to check that collapsed_pos + dr_zone_start >= DR_offset_map[*dr_iter] AND
						// that collapsed_pos < dr_zone_start - DR_offset_map[*dr_iter] + tmp_DR.length()
						//if(DR_offset_map[*dr_iter] <= dr_zone_start && dr_zone_start < (DR_offset_map[*dr_iter] + (int)(tmp_DR.length())) && collapsed_pos < (int)(tmp_DR.length()))
						if((collapsedPos + drAligner.getDRZoneStart() >= drAligner.offset(*dr_iter)) && (collapsedPos + drAligner.getDRZoneStart() - drAligner.offset(*dr_iter) < ((int)tmp_DR.length())))
						{
							// this is easy, we can compare based on this char only
							char decision_char = tmp_DR[drAligner.getDRZoneStart() - drAligner.offset(*dr_iter) + collapsedPos];
							collapsed_options2[decision_char] = collapsedOptions[decision_char];
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
					if(collapsedPos == 0)
					{
						logInfo("   ...ignoring (RLO SS)", 5);
					}
					else if(collapsedPos + drAligner.getDRZoneStart() == drAligner.getDRZoneEnd())
					{
						logInfo("   ...ignoring (RLO EE)", 5);
					}
					else
					{
						logInfo("   ...ignoring (RLO KK)", 5);
					}
#endif
					true_dr += drAligner.consensusAt(i);
					refinedDREnds[i] = true;
					collapsedOptions.clear();
				}
				else
				{
					// If it aint in collapsed_options2 it aint super fantastic!
					collapsedOptions.clear();
					collapsedOptions = collapsed_options2;
					
					// make the collapsed pos array specific and exit this loop
					collapsedPos += drAligner.getDRZoneStart();
					i = drAligner.getDRZoneEnd() + 1;
				}
			}
		}
	}
	logInfo("Consensus DR: " << true_dr, 1);
	return true_dr;
}

bool WorkHorse::parseGroupedDRs(int GID, int * nextFreeGID) 
{
	
    //-----
    // Cluster refinement and possible splitting for a Group ID
    //
    logInfo("Parsing group: " << GID, 4);
    		
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Find a Master DR for this group of DRs
    StringToken master_DR_token = -1;
    std::string master_DR_sequence = "**unset**";
    if(!findMasterDR(GID, &master_DR_token, &master_DR_sequence)) { return false; }
    
    
    // now we have the n most abundant kmers and one DR which contains them all
    // time to rock and rrrroll!
    
    Aligner dr_aligner((CRASS_DEF_CONS_ARRAY_RL_MULTIPLIER*mMaxReadLength), &mReads, &mStringCheck);
    dr_aligner.setMasterDR(master_DR_sequence);

    //++++++++++++++++++++++++++++++++++++++++++++++++
    // Set up the master DR's array and insert this guy into the main array
    populateCoverageArray(GID, dr_aligner );
    //++++++++++++++++++++++++++++++++++++++++++++++++
    // calculate consensus and diversity
	// use these variables to identify and store possible
	// collapsed clusters
	int collapsed_pos = -1;
	std::map<char, int> collapsed_options;            // holds the chars we need to split on
	std::map<int, bool> refined_DR_ends;              // so we can update DR ends based on consensus 
    std::string true_DR = calculateDRConsensus(GID, dr_aligner, *nextFreeGID, collapsed_pos, collapsed_options, refined_DR_ends);
    // check to make sure that the DR is not just some random long RE
    if((unsigned int)(true_DR.length()) > mOpts->highDRsize)
    {
        cleanGroup(GID);
        logInfo("Killed: {" << true_DR << "} cause' it was too long", 1);
        return false;
    }
    
    if (collapsed_options.size() == 0) {
        if((unsigned int)(true_DR.length()) < mOpts->lowDRsize)
        {
            cleanGroup(GID);
            logInfo("Killed: {" << true_DR << "} cause' the consensus was too short... (" << true_DR.length() << " ," << collapsed_options.size() << ")", 1);
            return false;
        }
        // QC the DR again for low complexity
        if (isRepeatLowComplexity(true_DR)) 
        {
            cleanGroup(GID);
            logInfo("Killed: {" << true_DR << "} cause' the consensus was low complexity...", 1);
            return false;
        }

        // test our true DR for highly abundant kmers
        try {
            float max_frequency;
            if (drHasHighlyAbundantKmers(true_DR, max_frequency) ) {
                cleanGroup(GID);
                logInfo("Killed: {" << true_DR << "} cause' the consensus contained highly abundant kmers: "<<max_frequency<<" > "<< CRASS_DEF_KMER_MAX_ABUNDANCE_CUTOFF, 1);
                return false;
            }
        } catch (crispr::exception& e) {
            std::cerr<<e.what()<<std::endl;
            throw crispr::runtime_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__, true_DR.c_str());
        }


        
        // update the DR start and ends
        // make a copy for ease of use
        int dr_zone_start = dr_aligner.getDRZoneStart(), dr_zone_end = dr_aligner.getDRZoneEnd();
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
        // set the updated values
        dr_aligner.setDRZoneStart(dr_zone_start);
        dr_aligner.setDRZoneEnd(dr_zone_end);
        
        char alphabet[4] = {'A', 'C', 'G', 'T'};

        // print out the consensus array
        if(isLogging(3))
        {
        	int show_xtra = 4;
        	int print_start = dr_aligner.getDRZoneStart() - show_xtra;
        	int print_end = dr_zone_end + show_xtra;
        	stringstream ss;
            ss << std::endl << "%, ";
            for(int i = print_start; i <= print_end; i++)
            {
            	if(i == dr_aligner.getDRZoneStart())
            		ss << "|, "; 
                ss << dr_aligner.conservationAt(i) << ", ";
            	if(i == dr_zone_end)
            		ss << "|, "; 
            }
            
            for(int j = 1; j <= 4; j++)
            {
                ss << std::endl;
                ss << alphabet[j-1] << ", ";
                for(int i = print_start; i <= print_end; i++)
                {
                	if(i == dr_aligner.getDRZoneStart())
                		ss << "|, "; 
                    ss << dr_aligner.coverageAt(i, alphabet[j-1]) << ", ";
                	if(i == dr_zone_end)
                		ss << "|, "; 
                }
            } 
            ss << std::endl << "$, ";
            for(int i = print_start; i <= print_end; i++)
            {
            	if(i == dr_aligner.getDRZoneStart())
            		ss << "|, "; 
                ss << dr_aligner.consensusAt(i) << ", ";
            	if(i == dr_zone_end)
            		ss << "|, "; 
            }
            logInfo(ss.str(), 3);
        }
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
            logInfo("Mapping \""<< co_iter->first << " : "  << co_iter->second << "\" to group: " << group << " "<< &mDR2GIDMap[group], 1);
            co_iter++;
        }
        
        DR_ClusterIterator dr_iter = (mDR2GIDMap[GID])->begin();
        while (dr_iter != (mDR2GIDMap[GID])->end()) 
        {
            std::string tmp_DR = mStringCheck.getString(*dr_iter);
            if(-1 != dr_aligner.offset(*dr_iter))
            {
                // check if the deciding character is within range of this DR
                if(dr_aligner.offset(*dr_iter) <= collapsed_pos && collapsed_pos < (dr_aligner.offset(*dr_iter) + (int)(tmp_DR.length())))
                {
                    logInfo("\tdeciding character within DR "<< *dr_iter,5);
                    // this is easy, we can compare based on this char only
                    char decision_char = tmp_DR[collapsed_pos - dr_aligner.offset(*dr_iter)];
                    (mDR2GIDMap[ coll_char_to_GID_map[ decision_char ] ])->push_back(*dr_iter);
                }
                else
                {
                    logInfo("\tRebuilding cluster from the ground up for DR " << *dr_iter,5);
                    // this is tricky, we need to completely break the group and re-cluster
                    // from the ground up based on reads
                    // get the offset from the start of the DR to the deciding char
                    // if it is negative, the dec char lies before the DR
                    // otherwise it lies after
                    int dec_diff = collapsed_pos - dr_aligner.offset(*dr_iter);
                    logInfo("\t\tBuilding form map",5);
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
                            logInfo("\t\tOne form found", 5);
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
                            logWarn("\t\tNo reads fit the form: " << tmp_DR, 1);
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
                            logInfo("\t\tMultiple forms, splitting reads", 5);
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
        logInfo("\tRemoving original group "<< GID, 5);
        cleanGroup(GID);
        
        logInfo("\tCalling the parser recursively", 4);
        
        // call this baby recursively with the new clusters
        std::map<char, int>::iterator cc_iter = coll_char_to_GID_map.begin();
        while(cc_iter != coll_char_to_GID_map.end())
        {
            parseGroupedDRs(cc_iter->second, nextFreeGID);
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
        logInfo("group: "<< GID<< " associated:" << &mDR2GIDMap[GID], 5);
        DR_ClusterIterator drc_iter = (mDR2GIDMap[GID])->begin();
        while(drc_iter != (mDR2GIDMap[GID])->end())
        {
        	logInfo("\tRepeat: "<<*drc_iter, 5);
            if(dr_aligner.offsetFind(*drc_iter) == dr_aligner.offsetEnd())
        	{
        		logError("1: Repeat "<< *drc_iter<<" in Group "<<GID <<" has no offset in DR_offset_map");
        	}
        	else
        	{
				if (dr_aligner.offset(*drc_iter) == -1 ) 
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
                        //if ((dr_aligner.offset(*drc_iter) - dr_aligner.getDRZoneStart()) < 0) {
                        //   // std::stringstream ss;
                        //    std::cerr << "front offset is a negative number\ndr_aligner.offset="<<dr_aligner.offset(*drc_iter)
                        //       << " dr_aligner.getDRZoneStart="<<dr_aligner.getDRZoneStart()
                        //       << " frontOffset="<<dr_aligner.offset(*drc_iter) - dr_aligner.getDRZoneStart();
                        //    //throw crispr::exception(__FILE__,
                        //    //                        __LINE__,
                        //    //                        __PRETTY_FUNCTION__,
                        //    //                        ss.str().c_str());
                        //} else {
                        //    std::cerr << "dr_aligner.offset="<<dr_aligner.offset(*drc_iter)
                        //       << " dr_aligner.getDRZoneStart="<<dr_aligner.getDRZoneStart()
                        //       << " frontOffset="<<dr_aligner.offset(*drc_iter) - dr_aligner.getDRZoneStart()<< std::endl;
                        //}
                        try {
                        //    std::cerr << "Alignment offset: "<< dr_aligner.offset(*drc_iter)<< " DR Zone Start: " <<dr_aligner.getDRZoneStart()<<std::endl;
						    (*read_iter)->updateStartStops((dr_aligner.offset(*drc_iter) - dr_aligner.getDRZoneStart()), &true_DR, mOpts);
                        } catch (crispr::exception &e) {
                            std::cerr <<dr_aligner.offset(*drc_iter) << " : "<<  dr_aligner.getDRZoneStart()<<std::endl;
                            logInfo("Dumping read set of group:", 1);
                            for (drc_iter = (mDR2GIDMap[GID])->begin(); drc_iter != (mDR2GIDMap[GID])->end(); drc_iter++) {
                                for (read_iter = mReads[*drc_iter]->begin(); read_iter != mReads[*drc_iter]->end(); read_iter++) {
                                    logInfoNoPrefix( *(*read_iter), 1); 
                                }
                            }
                            throw e;
                        }
						// reverse complement sequence if the true DR is not in its laurenized form
						if (rev_comp) 
						{
							try {
								(*read_iter)->reverseComplementSeq();
							} catch (crispr::exception& e) {
								std::cerr<<e.what()<<std::endl;
								throw crispr::exception(__FILE__,
								                        __LINE__,
								                        __PRETTY_FUNCTION__,
								                        "Failed to reverse complement sequence");
							}
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


void WorkHorse::cleanGroup(int GID)
{
    if(NULL != mDR2GIDMap[GID])
    {
        delete mDR2GIDMap[GID];
        mDR2GIDMap[GID] = NULL;
    }
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

bool WorkHorse::clusterDRReads(StringToken DRToken, 
                               int * nextFreeGID, 
                               std::map<std::string, int> * k2GIDMap, 
                               GroupKmerMap * groupKmerCountsMap)
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
	char ** kmers = NULL;
	int * kmer_offsets = NULL;
	try {
		kmers = new char*[num_mers];
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
			kmers[i] = new char [CRASS_DEF_KMER_SIZE+1];
		}
// use these offsets when we cut kmers, they are a component of the algorithm
		kmer_offsets = new int[num_mers];              
		for(int i = 0; i < num_mers; i++)
		{
			kmer_offsets[i] = i * -1; // Starts at [0, -1, -2, -3, -4, ...]
		}
	} catch(std::exception& e) {
		std::cerr<<"Attempting to alloc "<<CRASS_DEF_KMER_SIZE+1<<std::endl;
		throw crispr::exception(__FILE__, 
		                        __LINE__, 
		                        __PRETTY_FUNCTION__, 
		                        e.what());
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
    Vecstr homeless_kmers;
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
    Vecstr::iterator homeless_iter = homeless_kmers.begin();
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
#ifdef DEBUG
    //renderSpacerGraphs("Spacer_Preclean_");
#endif
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
#ifdef RENDERING
    std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering Debugging graphs using Graphviz"<<std::endl;
    logInfo("Rendering debug graphs" , 1);
#endif
    
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
                        std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering group "<<drg_iter->first<<std::endl;
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
#ifdef RENDERING
    std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering final spacer graphs using Graphviz"<<std::endl;
    logInfo("Rendering spacer graphs" , 1);
#endif
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
                NodeManager * current_manager = mDRs[mTrueDRs[drg_iter->first]];
                
                std::ofstream graph_file;

                
                std::string graph_file_prefix = mOpts->output_fastq + namePrefix + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first];
                std::string graph_file_name = graph_file_prefix + "_spacers.gv";
                
                // check to see if there is anything to print
                if ( current_manager->printSpacerGraph(graph_file_name, 
                                                       mTrueDRs[drg_iter->first], 
                                                       mOpts->longDescription, 
                                                       mOpts->showSingles))
                {
                    // add our group to the key
                    current_manager->printSpacerKey(key_file, 
                                                    10, 
                                                    namePrefix + to_string(drg_iter->first));
                    
                    // output the reads
                    std::string read_file_name = mOpts->output_fastq +  "Group_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first] + ".fa";

                    this->dumpReads(current_manager, read_file_name);
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
                    std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering group "<<drg_iter->first<<std::endl;
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

int WorkHorse::checkFileOrError(const char * fileName)
{
    // Test to see if the file is ok.
    struct stat inputDirStatus;
    int xStat = stat(fileName, &inputDirStatus);
    // stat failed
    if (xStat == -1) {
        return errno;
    }
    return 0;
}

bool WorkHorse::outputResults(std::string namePrefix)
{
	
    //-----
	// Print the cleaned? spacer graph, reads and the XML
	//

#ifdef RENDERING
    if(! mOpts->noRendering) {
        std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering final spacer graphs using Graphviz"<<std::endl;
        logInfo("Rendering spacer graphs" , 1);
    }
#endif
    // make a single file with all of the keys for the groups
    std::ofstream key_file;
    
    std::stringstream key_file_name;
    key_file_name << mOpts->output_fastq<<PACKAGE_NAME << "."<<mTimeStamp<<".keys.gv";
    key_file.open(key_file_name.str().c_str());
    
    if (!key_file) 
    {
        logError("Cannot open the key file");
        return 1;
    }
    
    gvGraphHeader(key_file, "Keys");

    
    
    // print all the assembly gossip to XML
	namePrefix += CRASS_DEF_CRISPR_EXT;
	logInfo("Writing XML output to \"" << namePrefix << "\"", 1);
	

    crispr::xml::writer * xml_doc = new crispr::xml::writer();
    int error_num;
    xercesc::DOMElement * root_element = xml_doc->createDOMDocument(CRASS_DEF_ROOT_ELEMENT, 
                                                                    CRASS_DEF_XML_VERSION, 
                                                                    error_num);
    
    if (!root_element && error_num) 
    {
        delete xml_doc;
        throw crispr::xml_exception(__FILE__,
                                    __LINE__,
                                    __PRETTY_FUNCTION__,
                                    "Unable to create xml document");
    }
    // go through the node managers and print the group info 
    // print all the inside information
    int final_out_number = 0;
    DR_Cluster_MapIterator drg_iter =  mDR2GIDMap.begin();
    for (drg_iter =  mDR2GIDMap.begin(); drg_iter != mDR2GIDMap.end(); drg_iter++) {

        // make sure that our cluster is real
        if (drg_iter->second == NULL) 
        {
            continue;
        }
        if(NULL == mDRs[mTrueDRs[drg_iter->first]])
        {
            continue;
        }

        NodeManager * current_manager = mDRs[mTrueDRs[drg_iter->first]];
        
        std::ofstream graph_file;        
        
        std::string graph_file_prefix = mOpts->output_fastq + "Spacers_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first];
        std::string graph_file_name = graph_file_prefix + "_spacers.gv";
        
        // check to see if there is anything to print
        if ( current_manager->printSpacerGraph(graph_file_name, 
                                               mTrueDRs[drg_iter->first], 
                                               mOpts->longDescription, 
                                               mOpts->showSingles))
        {
#ifdef RENDERING
            if (!mOpts->noRendering) 
            {
                // create a command string and call graphviz to make the image file
                std::cout<<"["<<PACKAGE_NAME<<"_imageRenderer]: Rendering group "<<drg_iter->first<<std::endl;
                std::string cmd = mOpts->layoutAlgorithm + " -Teps " + graph_file_name + " > "+ graph_file_prefix + ".eps";
                if(system(cmd.c_str()))
                {
                    logError("Problem running "<<mOpts->layoutAlgorithm<<" when rendering spacer graphs");
                    return 1;
                }
            }
#endif
            // add our group to the key
            current_manager->printSpacerKey(key_file, 
                                            10, 
                                            namePrefix + to_string(drg_iter->first));
            
            // output the reads
            std::string read_file_name = mOpts->output_fastq +  "Group_" + to_string(drg_iter->first) + "_" + mTrueDRs[drg_iter->first] + ".fa";
            this->dumpReads(current_manager, read_file_name, true);
            
            /* 
             *   Output the xml data to crass.crispr
             */
            std::string gid_as_string = "G" + to_string(drg_iter->first);
            final_out_number++;
            xercesc::DOMElement * group_elem = xml_doc->addGroup(gid_as_string, 
                                                                 mTrueDRs[drg_iter->first], 
                                                                 root_element);
            /*
             * <data> section
             */
            this->addDataToDOM(xml_doc, group_elem, drg_iter->first);
            
            /*
             * <metadata> section
             */
            this->addMetadataToDOM(xml_doc, group_elem, drg_iter->first);
            
            /*
             * <assembly> section
             */
            xercesc::DOMElement * assem_elem = xml_doc->addAssembly(group_elem);
            current_manager->printAssemblyToDOM(xml_doc, assem_elem, false);
        }
        else 
        {
            // should delete this guy since there are no spacers
            delete mDRs[mTrueDRs[drg_iter->first]];
            mDRs[mTrueDRs[drg_iter->first]] = NULL;
        }
    }
    std::cout<<"["<<PACKAGE_NAME<<"_graphBuilder]: "<<final_out_number<<" CRISPRs found!"<<std::endl;
    xml_doc->printDOMToFile(namePrefix);

    delete xml_doc;
    
    gvGraphFooter(key_file);
    key_file.close();
	return 0;
}

bool WorkHorse::addDataToDOM(crispr::xml::writer * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    try 
    {
        xercesc::DOMElement * data_elem = xmlDoc->addData(groupElement);
        if ((mDRs[mTrueDRs[groupNumber]])->haveAnyFlankers()) {
            xmlDoc->createFlankers(data_elem);
        }
        
        xercesc::DOMElement * sources_tag = data_elem->getFirstElementChild();
        std::set<StringToken> all_sources;
        
        for (xercesc::DOMElement * currentElement = data_elem->getFirstElementChild(); currentElement != NULL; currentElement = currentElement->getNextElementSibling()) 
        {
            if( xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->tag_Drs()))
            {
                // TODO: current implementation in Crass only supports a single DR for a group
                // in the future this will change, but for now ok to keep as a constant
                std::string drid = "DR1";
                xmlDoc->addDirectRepeat(drid, mTrueDRs[groupNumber], currentElement);
            }
            else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->tag_Spacers()))
            {
                // print out all the spacers for this group
                (mDRs[mTrueDRs[groupNumber]])->addSpacersToDOM(xmlDoc, currentElement, false, all_sources);
                
            }
            else if (xercesc::XMLString::equals(currentElement->getTagName(), xmlDoc->tag_Flankers()))
            {
                // should only get in here if there are flankers for the group

                // print out all the flankers for this group
                (mDRs[mTrueDRs[groupNumber]])->addFlankersToDOM(xmlDoc, currentElement, false, all_sources);

            }
        }
        (mDRs[mTrueDRs[groupNumber]])->generateAllsourceTags(xmlDoc, all_sources, sources_tag);
        
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

bool WorkHorse::addMetadataToDOM(crispr::xml::writer * xmlDoc, xercesc::DOMElement * groupElement, int groupNumber)
{
    try {
        
        std::stringstream notes;
        notes << "Run on "<< mTimeStamp;
        xercesc::DOMElement * metadata_elem = xmlDoc->addMetaData(groupElement);
        xercesc::DOMElement * prog_elem = xmlDoc->addProgram(metadata_elem);
        xmlDoc->addProgName(PACKAGE_NAME, prog_elem);
        xmlDoc->addProgVersion(PACKAGE_VERSION, prog_elem);
        xmlDoc->addProgCommand(mCommandLine, prog_elem);
        xmlDoc->addNotesToMetadata(notes.str(), metadata_elem);
        
        std::string file_name;
        char buf[4096];
        if(getcwd(buf, 4096) == NULL) {
        	crispr::exception(__FILE__,
                              __LINE__,
                              __PRETTY_FUNCTION__,
                              "Something went wrong getting the the CWD");
        }
        std::string absolute_dir = buf;
        absolute_dir += "/";
        // add in files if they exist
        if (!mOpts->logToScreen) 
        {
            // we whould have a log file
            file_name = mOpts->output_fastq + PACKAGE_NAME + "." + mTimeStamp + ".log";
            if (! checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("log", absolute_dir + file_name, metadata_elem);
            }
            else
            {
                throw crispr::no_file_exception(__FILE__, 
                                                __LINE__, 
                                                __PRETTY_FUNCTION__,
                                                (absolute_dir + file_name).c_str());
            }
        }
        
        
#ifdef DEBUG
        // check for debuging .gv files
        if (!mOpts->noDebugGraph) 
        {
            file_name = mOpts->output_fastq + "Group_"; 
            std::string file_sufix = to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + "_debug.gv";
            if (! checkFileOrError((file_name + file_sufix).c_str())) 
            {
                xmlDoc->addFileToMetadata("data", absolute_dir + file_name + file_sufix, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, 
                                                __LINE__, 
                                                __PRETTY_FUNCTION__, 
                                                (absolute_dir + file_name + file_sufix).c_str() );
            }
            
            // and now for the cleaned .gv
            file_name = mOpts->output_fastq + "Clean_";
            if (! checkFileOrError((file_name + file_sufix).c_str())) 
            {
                xmlDoc->addFileToMetadata("data", absolute_dir + file_name + file_sufix, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, 
                                                __LINE__, 
                                                __PRETTY_FUNCTION__, 
                                                (absolute_dir + file_name + file_sufix).c_str() );                
            }
        }

        
#endif // DEBUG
        
#ifdef RENDERING
        // check for image files
#ifdef DEBUG
        if (!mOpts->noDebugGraph) 
        {
            file_name = mOpts->output_fastq + "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
            if (! checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("image", absolute_dir + file_name, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
            
            file_name = mOpts->output_fastq + "Clean_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".eps";
            
            if (! checkFileOrError(file_name.c_str())) 
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
            if (! checkFileOrError(file_name.c_str())) 
            {
                xmlDoc->addFileToMetadata("image", absolute_dir + file_name, metadata_elem);
            } 
            else 
            {
                throw crispr::no_file_exception(__FILE__, __LINE__, __PRETTY_FUNCTION__,(absolute_dir + file_name).c_str());
            }
        } 


#endif // RENDERING
        
        // add in the final Spacer graph
        file_name = mOpts->output_fastq + "Spacers_"; 
        std::string file_sufix = to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + "_spacers.gv";
        if (! checkFileOrError((file_name + file_sufix).c_str())) 
        {
            xmlDoc->addFileToMetadata("data", absolute_dir + file_name + file_sufix, metadata_elem);
        } 
        else 
        {
            throw crispr::no_file_exception(__FILE__, 
                                            __LINE__, 
                                            __PRETTY_FUNCTION__, 
                                            (absolute_dir + file_name + file_sufix).c_str() );
        }

        
        // check the sequence file
        file_name = mOpts->output_fastq +  "Group_" + to_string(groupNumber) + "_" + mTrueDRs[groupNumber] + ".fa";
        if (! checkFileOrError(file_name.c_str())) 
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



