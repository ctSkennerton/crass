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
    typedef std::map<int, std::map<std::string, int> * > GroupKmerMap;

    class Storage {
    private:
        std::vector<crass::RawRead>                             mReads;
        std::unordered_map<StringToken, std::list<size_t> >     mRepeatsToReads;
        std::unordered_map<int, std::list<StringToken> >        mRepeatCluster;
        std::list<StringToken>                                  mNonRedundantRepeats;
        int                                                     mNextFreeGID;
        StringCheck                                             mRepeatTokenizer;
        
        void linkReadToRepeat(const size_t readPosition, const StringToken repeat);
        void linkReadToRepeat(const size_t read1Position, const size_t read2Position, const StringToken repeat);
        StringToken tokenizeRepeat( std::string& repeat);
        bool orientateLowLexi(crass::RawRead& read);
        void removeRedundantRepeats(std::vector<std::string>& repeatVector);
        void createNonRedundantSet(GroupKmerMap& groupKmerCountsMap);
    public:
        Storage(){mNextFreeGID = 0;}
        
        void add(crass::RawRead& read);
        void add(crass::RawRead& read1, crass::RawRead& read2);
        void clusterRepeats(int minKmerCount);
        int clusterRepeats(std::string& outputDirectory, float identityThreshold, int threads);
        void inspect(std::ostream& out);
    };
    
}
#endif /* defined(__crass__Storage__) */
