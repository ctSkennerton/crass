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
using namespace std;
class PatternMatcher{
public:
    static int bmpSearch(const string& text, const string& pattern);
    
private:
    static vector<int> computeBmpLast(const string& pattern);
    
    PatternMatcher();
    PatternMatcher(const PatternMatcher&);
    const PatternMatcher& operator=(const PatternMatcher&);
};
#endif 