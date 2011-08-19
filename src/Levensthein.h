/*
 *  Levensthein.h is part of the MCD project
 *  
 *  This code was downloaded from http://www.merriampark.com/ldcpp.htm
 *  Copyright (c) Merriam Park Software 2009
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

#ifndef __LEV_H
#define __LEV_H
#include <string>
#include <vector>

typedef std::vector< std::vector<int> > Tmatrix; 


int LevenstheinDistance( std::string source,  std::string target);
#endif
