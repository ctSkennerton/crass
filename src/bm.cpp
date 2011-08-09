/*
 *  bm.cpp is part of the CRisprASSembler project
 *  
 *  This code was downloaded from http://dev-faqs.blogspot.com/2010/05/boyer-moore-algorithm.html
 *  Copyright (c) 2010 dev-faqs.blogspot.com
 *
 *                     A B R A K A D A B R A
 *                      A B R A K A D A B R
 *                       A B R A K A D A B
 *                        A B R A K A D A       	
 *                         A B R A K A D
 *                          A B R A K A
 *                           A B R A K
 *                            A B R A
 *                             A B R
 *                              A B
 *                               A
 */


#include "bm.h"
#include <algorithm>

using namespace std;

int PatternMatcher::bmpSearch(const std::string &text, const std::string &pattern){
    size_t textSize = text.size();
    size_t patternSize = pattern.size();
    if(textSize == 0 || patternSize == 0){
        return -1;
    }
    if(patternSize > textSize){
        return -1;
    }
    
    vector<int> bmpLast = computeBmpLast(pattern);
    size_t tIdx = patternSize - 1;
    size_t pIdx = patternSize - 1;
    while(tIdx < textSize)
    {
        if(pattern[pIdx] == text[tIdx])
        {
            if(pIdx == 0)
            {   //found a match
                return (int)tIdx;
            }
            tIdx--;
            pIdx--;
        }
        else 
        {
            //Character Jump Heuristics
            int lastOccur = bmpLast[text[tIdx]];
            tIdx = tIdx + patternSize - min<int>((int)pIdx, 1 + lastOccur);
            pIdx = patternSize - 1;
        }
    }
    return - 1;
}

// slight variation of the code above so that when a match is found it is pushed to 
// a vector of starting positions
void PatternMatcher::bmpMultiSearch(const std::string &text, const std::string &pattern, std::vector<int> &startOffsetVec )
{
    size_t textSize = text.size();
    size_t patternSize = pattern.size();
    if(textSize == 0 || patternSize == 0){
        return;
    }
    if(patternSize > textSize){
        return;
    }
    
    vector<int> bmpLast = computeBmpLast(pattern);
    size_t tIdx = patternSize - 1;
    size_t pIdx = patternSize - 1;
    while(tIdx < textSize)
    {
        if(pattern[pIdx] == text[tIdx])
        {
            if(pIdx == 0)
            {   //found a match
                startOffsetVec.push_back((int)tIdx) ;
            }
            tIdx--;
            pIdx--;
        }
        else 
        {
            //Character Jump Heuristics
            int lastOccur = bmpLast[text[tIdx]];
            tIdx = tIdx + patternSize - min<int>((int)pIdx, 1 + lastOccur);
            pIdx = patternSize - 1;
        }
    }
}


vector<int> PatternMatcher::computeBmpLast(const std::string &pattern){
    const size_t NUM_ASCII_CHARS = 128;
    vector<int> bmpLast(NUM_ASCII_CHARS);
    for(size_t i = 0; i < NUM_ASCII_CHARS; i++){
        bmpLast[i] = -1;
    }
    for(size_t i = 0; i < pattern.size(); i++){
        bmpLast[pattern[i]] = (int)i;
    }
    return bmpLast;
}
