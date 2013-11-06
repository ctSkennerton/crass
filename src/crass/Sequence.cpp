//
//  Sequence.cpp
//  crass
//
//  Created by Connor Skennerton on 1/11/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#include "Sequence.h"
#include <algorithm>

using namespace crass;


char comp_tab[] = {
    0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
    16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
    32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
    48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
    64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

std::string crass::reverseComplement(std::string str)
{
    
	int l = static_cast<int>(str.length());
    char * revcomp_str = new char[l+1];
    for (int i = 0; i <=l; i++) {
		revcomp_str[i]='\0';
    }
    int i, c0, c1;
    for (i = 0; i < l>>1; ++i)
    {
        c0 = comp_tab[(int)str[i]];
        c1 = comp_tab[(int)str[l - 1 - i]];
        revcomp_str[i] = c1;
        revcomp_str[l - 1 - i] = c0;
    }
    if (l&1)
    {
        revcomp_str[l>>1] = comp_tab[(int)str[l>>1]];
    }
    
    std::string ret = revcomp_str;
    delete [] revcomp_str;
    
    return ret;
}

std::string crass::laurenize (std::string seq1)
{
    std::string seq2 = reverseComplement(seq1);
    if (seq1 < seq2)
    {
        return seq1;
    }
    return seq2;
}

void RawRead::revComp() {
    mSeq = reverseComplement(mSeq);
    mRepeatPositions.reverseRepeatPositions(mSeq.size() - 1);
    if(!mQual.empty()) {
        std::reverse(mQual.begin(), mQual.end());
    }
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
    switch (numberOfRepeats()) {
        case 0: {
            // error
            return -1;
            break;
        }
        case 1: {
            return mRepeatPositions.repeatLengthAt(0);
            break;
        }
        case 2: {
            int first_length = mRepeatPositions.repeatLengthAt(0);
            int second_length = mRepeatPositions.repeatLengthAt(1);
            bool sp = startPartial();
            bool ep = endPartial();
            if (sp && ep) {
                return (first_length > second_length) ? first_length : second_length;
            } else if (sp) {
                return second_length;
            } else {
                return first_length;
            }
            break;
        }
        default: {
            auto it = repeatBegin();
            if(startPartial()) {
                return mRepeatPositions.repeatLengthAt(1);
            } else {
                return mRepeatPositions.repeatLengthAt(2);
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