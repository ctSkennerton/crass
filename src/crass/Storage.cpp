//
//  Storage.cpp
//  crass
//
//  Created by Connor Skennerton on 20/10/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#include "Storage.h"
#include "SeqUtils.h"

using namespace crass;

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

void Storage::linkReadToRepeat(const size_t readPosition, const StringToken repeat){
    mRepeatsToReads[repeat].push_back(readPosition);
}

void Storage::linkReadToRepeat(const size_t read1Position, const size_t read2Position, const StringToken repeat){
    mRepeatsToReads[repeat].push_back(read1Position);
    mRepeatsToReads[repeat].push_back(read2Position);
}

StringToken Storage::tokenizeRepeat( std::string &repeat) {
    StringToken t = mRepeatTokenizer.getToken(repeat);
    if(t == 0) {
        return mRepeatTokenizer.addString(repeat);
    } else {
        return t;
    }
}

bool Storage::orientateLowLexi(RawRead& read)
{
    std::string tmp_dr = *(read.getFirstNonPartialRepeat());
    std::string rev_comp = reverseComplement(tmp_dr);        
    if (tmp_dr < rev_comp)
    {
        return true;
    }
    else
    {
        read.revComp();
        return false;
    }
}


void Storage::add(RawRead& read) {
    orientateLowLexi(read);
    mReads.push_back(read);
    
    std::string s = *(read.getFirstNonPartialRepeat());
    StringToken t = tokenizeRepeat(s);
    linkReadToRepeat(mReads.size() - 1, t);
}

void Storage::add(RawRead& read1, RawRead& read2) {
    
    // assume that both reads enter in the same orientation
    // i.e. FF or RR
    // do not allow for FR or RF.  this shouldn't happen since
    // the search algorithms need them in the same orientation.
    
    bool was_reversed = orientateLowLexi(read1);
    if (was_reversed) {
        read2.revComp();
    }
    
    mReads.push_back(read1);
    mReads.push_back(read2);

    auto read2_it = mReads.end() - 1;
    auto read1_it = read2_it - 1;

    (*read2_it).previousRead(&(*read1_it));
    
    (*read1_it).nextRead(&(*read2_it));
    
    std::string s = *(read1.getFirstNonPartialRepeat());
    StringToken t = tokenizeRepeat(s);
    linkReadToRepeat(mReads.size() - 2, mReads.size() - 1, t);
    
}


void Storage::clusterRepeats(int minKmerCount)
{
    logInfo("Reducing list of potential DRs (1): Initial clustering", 1);
    logInfo("Reticulating splines...", 1);
    
    const int kmer_size = 11;
    
    std::map<std::string, int> k2GIDMap;
    GroupKmerMap groupKmerCountsMap;
    
    for( auto it = mRepeatsToReads.begin(); it != mRepeatsToReads.end(); ++it) {
        std::string repeat = mRepeatTokenizer.getString(it->first);
        int num_mers;
        std::vector<std::string> kmers;
        crass::cutIntoKmers(repeat, kmer_size, kmers);

        std::vector<std::string> homeless_kmers;
        std::map<int, int> group_count;
        std::map<std::string, int> local_kmer_CountMap;
        
        int group = 0;
        for(auto it2 = kmers.begin() ; it2 != kmers.end(); ++it2)
        {
            // make it a string!
            
            std::string tmp_str = laurenize(*it2);
            
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
            std::map<std::string, int>::iterator k2g_iter = k2GIDMap.find(tmp_str);
            if(k2g_iter == k2GIDMap.end())
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
                        if(minKmerCount <= group_count[k2g_iter->second])
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
            group = mNextFreeGID++;
            groupKmerCountsMap[group] = new std::map<std::string, int>;
        }
        
        // we need to record the group for this mofo!
        mRepeatCluster[group].push_back(it->first);
        
        // we need to assign all homeless kmers to the group!
        std::vector<std::string>::iterator homeless_iter = homeless_kmers.begin();
        while(homeless_iter != homeless_kmers.end())
        {
            k2GIDMap[*homeless_iter] = group;
            homeless_iter++;
        }
        
        // we need to fix up the group counts
        std::map<std::string, int>::iterator local_count_iter = local_kmer_CountMap.begin();
        while(local_count_iter != local_kmer_CountMap.end())
        {
            (*groupKmerCountsMap[group])[local_count_iter->first] += local_count_iter->second;
            local_count_iter++;
        }

    }
    
    createNonRedundantSet(groupKmerCountsMap);
    
    GroupKmerMap::iterator group_count_iter;
    for(group_count_iter =  groupKmerCountsMap.begin();
        group_count_iter != groupKmerCountsMap.end();
        ++group_count_iter)
    {
        // delete the kmer count lists cause we're finsihed with them now
        if(NULL != group_count_iter->second)
        {
            delete group_count_iter->second;
            group_count_iter->second = NULL;
        }
    }
        
}

void Storage::removeRedundantRepeats(std::vector<std::string>& repeatVector)
{
    // given a vector of repeat sequences, will order the vector based on repeat
    // length and then remove longer repeats if there is a shorter one that is
    // a perfect substring
    std::sort(repeatVector.begin(), repeatVector.end(), sortLengthAssending);
    std::vector<std::string>::iterator iter;
    
    // go though all of the patterns and determine which are substrings
    // clear the string if it is
    for (iter = repeatVector.begin(); iter != repeatVector.end(); iter++) {
        if (iter->empty()) {
            continue;
        }
        std::vector<std::string>::iterator iter2;
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
    std::vector<std::string>::iterator empty_iter = std::partition(repeatVector.begin(), repeatVector.end(), isNotEmpty);
    // remove all the empties from the list
    repeatVector.erase(empty_iter, repeatVector.end());
}


void Storage::createNonRedundantSet(GroupKmerMap& groupKmerCountsMap)
{

    std::cout<<'['<<PACKAGE_NAME<<"_clusterCore]: "<<mRepeatsToReads.size()<<" variants mapped to "<<mRepeatCluster.size()<<" clusters"<<std::endl;
    std::cout<<'['<<PACKAGE_NAME<<"_clusterCore]: creating non-redundant set"<<std::endl;

    for(auto cluster_it = mRepeatCluster.begin(); cluster_it != mRepeatCluster.end(); ++cluster_it)
    {
        logInfo("-------------", 4);
        logInfo("Group: " << cluster_it->first, 4);
        std::vector<std::string> clustered_repeats;
        for(auto repeat_it = (cluster_it->second).begin(); repeat_it != (cluster_it->second).end(); ++repeat_it) {
            std::string tmp = mRepeatTokenizer.getString(*repeat_it);
            clustered_repeats.push_back(tmp);
            logInfo(tmp, 4);
        }
        
        logInfo("-------------", 4);
        removeRedundantRepeats(clustered_repeats);
        std::vector<std::string>::iterator cr_iter;
        for (cr_iter = clustered_repeats.begin(); cr_iter != clustered_repeats.end(); ++cr_iter) {
            mNonRedundantRepeats.push_back(mRepeatTokenizer.getToken(*cr_iter));
        }
    }
}


void Storage::inspect(std::ostream &out) {
    out <<"Number of reads: "<<mReads.size()<<std::endl;
    out << "Number of patterns: "<<mRepeatsToReads.size()<<std::endl;
    out << "Number of clusters: "<<mRepeatCluster.size()<<std::endl;

    if(!mRepeatCluster.empty()) {
        for (auto it = mRepeatCluster.begin(); it != mRepeatCluster.end(); ++it) {
            out<< it->first<<"\t"<<(it->second).size()<<std::endl;
        }
    }
}