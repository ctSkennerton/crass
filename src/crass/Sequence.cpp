//
//  Sequence.cpp
//  crass
//
//  Created by Connor Skennerton on 1/11/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#include "Sequence.h"
#include "LoggerSimp.h"
#include <algorithm>

using namespace crass;


void RawRead::revComp() {
    mSeq = reverseComplement(mSeq);
    mRepeatPositions.reverseRepeatPositions(mSeq.size() - 1);
    if(!mQual.empty()) {
        std::reverse(mQual.begin(), mQual.end());
    }
    mRevComp = !mRevComp;
    mRepeatLowLexi = !mRepeatLowLexi;
}
    
std::string RawRead::orientateRepeatLowLexi(void)
{
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    int num_repeats = mRepeatPositions.numberOfRepeats();
    if (num_repeats == 1)
    {
        
        tmp_dr = *(repeatStringBegin());
        rev_comp = reverseComplement(tmp_dr);
    }
    else if (2 == num_repeats)
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // take the second
        RepeatArray::RepeatIterator it = mRepeatPositions.repeatAt(0);
        RepeatArray::RepeatIterator it1 = mRepeatPositions.repeatAt(1);

        if ((*it).first == 0)
        {
            tmp_dr = (*repeatStringAt(1));
            rev_comp = reverseComplement(tmp_dr);
        }
        else if ((*it1).second == static_cast<int>(mSeq.length()))
        {
            tmp_dr = (*repeatStringAt(0));
            rev_comp = reverseComplement(tmp_dr);

        }
        else
        {
            // if they both are then just take whichever is longer
            int lenA = (*it).second - (*it).first;
            int lenB = (*it1).second - (*it1).first;
            
            if (lenA > lenB)
            {
                tmp_dr = (*repeatStringAt(0));
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = (*repeatStringAt(1));
                rev_comp = reverseComplement(tmp_dr);
            }
        }
    }
    // long read more than two repeats
    else
    {
        // take the second
        tmp_dr = (*repeatStringAt(1));
        rev_comp = reverseComplement(tmp_dr);
        
    }
    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        mRepeatLowLexi = true;
#ifdef DEBUG
        //logInfo("DR in low lexi"<<mSeq, 9);
#endif
        return tmp_dr;
    }
    else
    {
        revComp();
        mRepeatLowLexi = false;
#ifdef DEBUG
        //logInfo("DR not in low lexi"<<mSeq, 9);
#endif
        return rev_comp;
    }
}

int RawRead::getFirstNonPartialRepeatLength() {
    RepeatStringIterator it = getFirstNonPartialRepeat();
    return static_cast<int>((*it).length());
}

RawRead::RepeatStringIterator crass::RawRead::getFirstNonPartialRepeat() {
    switch (numberOfRepeats()) {
        case 0: {
            // error
            return repeatStringAt(0);
            break;
        }
        case 1: {
            return repeatStringAt(0);
            break;
        }
        case 2: {
            int first_length = mRepeatPositions.repeatLengthAt(0);
            int second_length = mRepeatPositions.repeatLengthAt(1);
            bool sp = startPartial();
            bool ep = endPartial();
            if (sp && ep) {
                return (first_length > second_length) ? repeatStringAt(0) : repeatStringAt(1);
            } else if (sp) {
                return repeatStringAt(1);
            } else {
                return repeatStringAt(0);
            }
            break;
        }
        default: {
            auto it = repeatBegin();
            if(startPartial()) {
                return repeatStringAt(1);
            } else {
                return repeatStringAt(0);
            }
            break;
        }
    }
}

RepeatArray::RepeatIterator crass::RawRead::getFirstNonPartialRepeatPositions() {
    switch (numberOfRepeats()) {
        case 0: {
            // error
            return repeatAt(0);
            break;
        }
        case 1: {
            return repeatAt(0);
            break;
        }
        case 2: {
            int first_length = mRepeatPositions.repeatLengthAt(0);
            int second_length = mRepeatPositions.repeatLengthAt(1);
            bool sp = startPartial();
            bool ep = endPartial();
            if (sp && ep) {
                return (first_length > second_length) ? repeatAt(0) : repeatAt(1);
            } else if (sp) {
                return repeatAt(1);
            } else {
                return repeatAt(0);
            }
            break;
        }
        default: {
            auto it = repeatBegin();
            if(startPartial()) {
                return repeatAt(1);
            } else {
                return repeatAt(0);
            }
            break;
        }
    }
}

#ifdef crass_RawRead_main
int main() {
    
    RepeatArray r = RepeatArray();
    RawRead seq1 = RawRead("1", "my seq", "AAAAAAAGGGGGGG", "@@#$%3&^%#@ABDH", r);
    seq1.push_back(0, 3);
    seq1.push_back(7, 12);
    seq1.inspect();
    std::cout<<"Testing reverse complement:"<<std::endl;
    seq1.revComp();
    seq1.inspect();
    seq1.revComp();
    std::cout<<"Testing DR LowLexi orientation:"<<std::endl;
    std::string s = seq1.orientateRepeatLowLexi();
    std::cout<<s<<std::endl;
    seq1.inspect();
    std::cout<<"Testing repeat string iteration:"<<std::endl;
    for(auto it = seq1.repeatStringBegin(); it != seq1.repeatStringEnd(); ++it) {
        std::cout<<*it<<std::endl;
    }
    std::cout<<"Testing spacer string iteration:"<<std::endl;
    for(auto it = seq1.spacerStringBegin(); it != seq1.spacerStringEnd(); ++it) {
        std::cout<<*it<<std::endl;
    }
    return 0;
}
#endif