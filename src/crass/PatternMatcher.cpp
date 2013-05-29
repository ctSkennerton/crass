/*
 *  PatternMatcher.cpp is part of the CRisprASSembler project
 *  Levensthein code was downloaded from http://www.merriampark.com/ldcpp.htm
 *  Copyright (c) Merriam Park Software 2009
 *  
 *  Boyer-Moore code was downloaded from http://dev-faqs.blogspot.com/2010/05/boyer-moore-algorithm.html
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


#include "PatternMatcher.h"
#include <algorithm>

int PatternMatcher::bmpSearch(const std::string &text, const std::string &pattern){
    size_t textSize = text.size();
    size_t patternSize = pattern.size();
    if(textSize == 0 || patternSize == 0){
        return -1;
    }
    if(patternSize > textSize){
        return -1;
    }
    
    std::vector<int> bmpLast = computeBmpLast(pattern);
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
            tIdx = tIdx + patternSize - std::min<int>((int)pIdx, 1 + lastOccur);
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
    
    std::vector<int> bmpLast = computeBmpLast(pattern);
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
            tIdx = tIdx + patternSize - std::min<int>((int)pIdx, 1 + lastOccur);
            pIdx = patternSize - 1;
        }
    }
}


std::vector<int> PatternMatcher::computeBmpLast(const std::string &pattern){
    const size_t NUM_ASCII_CHARS = 128;
    std::vector<int> bmpLast(NUM_ASCII_CHARS);
    for(size_t i = 0; i < NUM_ASCII_CHARS; i++){
        bmpLast[i] = -1;
    }
    for(size_t i = 0; i < pattern.size(); i++){
        bmpLast[pattern[i]] = (int)i;
    }
    return bmpLast;
}

int PatternMatcher::levenstheinDistance( std::string& source,  std::string& target) {
    
    // Step 1
    
    int n = (int)source.length();
    int m = (int)target.length();
    if (n == 0) {
        return m;
    }
    if (m == 0) {
        return n;
    }
    
    // Good form to declare a TYPEDEF
    
    
    Tmatrix matrix(n+1);
    
    // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
    // allow for allocation on declaration of 2.nd dimension of vec of vec
    
    for (int i = 0; i <= n; i++) {
        matrix[i].resize(m+1);
    }
    
    // Step 2
    
    for (int i = 0; i <= n; i++) {
        matrix[i][0]=i;
    }
    
    for (int j = 0; j <= m; j++) {
        matrix[0][j]=j;
    }
    
    // Step 3
    
    for (int i = 1; i <= n; i++) {
        
        char s_i = source[i-1];
        
        // Step 4
        
        for (int j = 1; j <= m; j++) {
            
            char t_j = target[j-1];
            
            // Step 5
            
            int cost;
            if (s_i == t_j) {
                cost = 0;
            }
            else {
                cost = 1;
            }
            
            // Step 6
            
            int above = matrix[i-1][j];
            int left = matrix[i][j-1];
            int diag = matrix[i-1][j-1];
            int cell = std::min( above + 1, std::min(left + 1, diag + cost));
            
            // Step 6A: Cover transposition, in addition to deletion,
            // insertion and substitution. This step is taken from:
            // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's 
            // Enhanced Dynamic Programming ASM Algorithm"
            // (http://www.acm.org/~hlb/publications/asm/asm.html)
            
            if (i>2 && j>2) {
                int trans=matrix[i-2][j-2]+1;
                if (source[i-2]!=t_j) trans++;
                if (s_i!=target[j-2]) trans++;
                if (cell>trans) cell=trans;
            }
            
            matrix[i][j]=cell;
        }
    }
    
    // Step 7
    
    return matrix[n][m];
}

float PatternMatcher::getStringSimilarity(std::string& s1, std::string& s2)
{
    float max_length = std::max(s1.length(), s2.length());
    if(/*max_length > 10 ||*/ s1.length() < 3 || s2.length() < 3)
    	return 0;
    float edit_distance =  levenstheinDistance(s1 ,  s2);
    return 1.0 - (edit_distance/max_length);
}

