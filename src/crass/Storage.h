//
//  Storage.h
//  crass
//
//  Created by Connor Skennerton on 20/10/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#ifndef __crass__Storage__
#define __crass__Storage__

#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <map>

#include "Sequence.h"
#include "StringCheck.h"
#include "LoggerSimp.h"


namespace crass {
    
    class Search;
    class ReadIterator;

    //typedef std::map<int, std::map<std::string, int> * > GroupKmerMap;
    
    // key is an integer representing a group
    // value in a list of StringToken that are the clustered repeats
    typedef std::unordered_map<int, std::list<StringToken> > RepeatCluster_t;
    
    // key is a StringToken that represents a repeat sequence
    // value is a list of indexes into the vector containing all of the RawRead objects
    typedef std::unordered_map<StringToken, std::list<size_t> > ReadCluster_t;

    class Storage {
        friend class Search;
        
    protected:
        std::vector<crass::RawRead>                             mReads;
        ReadCluster_t                                           mRepeatsToReads;
        RepeatCluster_t                                         mRepeatCluster;
        std::list<StringToken>                                  mNonRedundantRepeats;
        int                                                     mNextFreeGID;
        StringCheck                                             mRepeatTokenizer;
        
        void linkReadToRepeat(const size_t readPosition, const StringToken repeat);
        void linkReadToRepeat(const size_t read1Position, const size_t read2Position, const StringToken repeat);
        StringToken tokenizeRepeat( std::string& repeat);
        bool orientateLowLexi(crass::RawRead& read);
        void removeRedundantRepeats(std::vector<std::string>& repeatVector);
        void createNonRedundantSet();
    public:
        /*
         
         */
        class ReadFromRepeatIterator {
        public:
            typedef ReadFromRepeatIterator self_t;
            typedef crass::RawRead value_t;
            typedef value_t& reference_t;
            typedef value_t* pointer_t;
            
            ReadFromRepeatIterator(std::vector<crass::RawRead>& reads, std::list<size_t>::iterator si) : data(reads), current_pos(si){}
            
            reference_t operator *() {
                return data.at(*current_pos);
            }
            pointer_t operator ->() {
                return &data.at(*current_pos);
            }
            
            self_t operator++() { self_t i = *this; ++current_pos; return i; }
            self_t operator++(int junk) { ++current_pos; return *this; }
            
            
            bool operator==( const self_t& rhs) { return current_pos == rhs.current_pos; }
            bool operator!=( const self_t& rhs) { return current_pos != rhs.current_pos; }

                
                
            private:
                std::vector<crass::RawRead>& data;
                std::list<size_t>::iterator current_pos;
            };
        
        Storage(){mNextFreeGID = 0;}
        
        void add(crass::RawRead& read);
        void add(crass::RawRead& read1, crass::RawRead& read2);
        void clusterRepeats(int minKmerCount);
        int clusterRepeats(std::string& outputDirectory, float identityThreshold, int threads);
        void inspect(std::ostream& out);
        void dumpReads(std::ostream& out);
        
        crass::RawRead& getRead(size_t i) {
            return mReads.at(i);
        }
        
        size_t numberOfReads(){return mReads.size();}
        std::vector<crass::RawRead>::iterator readBegin(){return mReads.begin();}
        std::vector<crass::RawRead>::iterator readEnd(){return mReads.end();}
        
        RepeatCluster_t::iterator repeatClusterBegin(){return mRepeatCluster.begin();}
        RepeatCluster_t::iterator repeatClusterEnd(){return mRepeatCluster.end();}
        
        ReadCluster_t::iterator readClusterBegin(){return mRepeatsToReads.begin();}
        ReadCluster_t::iterator readClusterEnd(){return mRepeatsToReads.end();}
        ReadFromRepeatIterator readFromRepeatBegin(StringToken repeat) { return ReadFromRepeatIterator(mReads, mRepeatsToReads[repeat].begin());}
        ReadFromRepeatIterator readFromRepeatEnd(StringToken repeat) { return ReadFromRepeatIterator(mReads, mRepeatsToReads[repeat].end());}

            
        void getNonRedundantRepeats(std::vector<std::string>& repeats);

    };
    
}
#endif /* defined(__crass__Storage__) */
