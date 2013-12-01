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

#include "Sequence.h"
#include "StringCheck.h"
#include "LoggerSimp.h"
namespace crass {
    
    class Storage {
    private:
        std::vector<crass::RawRead>           mReads;
        std::unordered_map<StringToken, std::list<size_t> >     mRepeatsToReads;
        StringCheck                           mRepeatTokenizer;
        
        void linkReadToRepeat(const size_t readPosition, const StringToken repeat);
        void linkReadToRepeat(const size_t read1Position, const size_t read2Position, const StringToken repeat);
        StringToken tokenizeRepeat( std::string& repeat);
        bool orientateLowLexi(crass::RawRead& read);

    public:
        void add(crass::RawRead& read);
        void add(crass::RawRead& read1, crass::RawRead& read2);
        void inspect(std::ostream& out);
    };
    
}
#endif /* defined(__crass__Storage__) */
