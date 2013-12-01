//
//  Storage.cpp
//  crass
//
//  Created by Connor Skennerton on 20/10/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#include "Storage.h"

using namespace crass;

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

void Storage::inspect(std::ostream &out) {
    out <<"Number of reads: "<<mReads.size()<<std::endl;
    out << "Number of patterns:"<<mRepeatsToReads.size()<<std::endl;
}