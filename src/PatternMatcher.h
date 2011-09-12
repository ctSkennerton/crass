/*
 *  bm.h is part of the CRisprASSembler project
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

#ifndef _BM_H_
#define _BM_H_
#include <iostream>
#include <string>
#include <vector>
typedef std::vector< std::vector<int> > Tmatrix; 

class PatternMatcher{
public:
    static int bmpSearch(const std::string& text, const std::string& pattern);
    
    static void bmpMultiSearch(const std::string &text, const std::string &pattern, std::vector<int> &startOffsetVec);
    
    static int levenstheinDistance( std::string& source,  std::string& target);
    
    static float getStringSimilarity(std::string& s1, std::string& s2);

private:
    static std::vector<int> computeBmpLast(const std::string& pattern);
    
    PatternMatcher();
    PatternMatcher(const PatternMatcher&);
    const PatternMatcher& operator=(const PatternMatcher&);
};
#endif 
