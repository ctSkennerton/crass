/*
 *  SeqUtils.cpp is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011 Connor Skennerton & Michael Imelfort. All rights reserved. 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
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


#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <config.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "SeqUtils.h"

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

std::string reverseComplement(std::string str)
{
    int l = (int)str.length();
    char * revcomp_str = new char[l+1];
    for (int i = 0; i <=l; i++) {
        revcomp_str[i]=NULL;
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

std::string laurenize (std::string seq1)
{
    std::string seq2 = reverseComplement(seq1);
    if (seq1 < seq2)
    {
        return seq1;
    }
    return seq2;
}


gzFile getFileHandle(const char * inputFile)
{
    gzFile fp;
    if ( strcmp(inputFile, "-") == 0 ) 
    {
        fp = gzdopen(fileno(stdin), "r");
    }
    else 
    {
        fp = gzopen(inputFile, "r");
    }
    
    if ( (fp == NULL) && (strcmp(inputFile, "-") != 0) ) 
    {
        std::cerr<< PACKAGE_NAME<<" : [ERROR] Could not open FASTQ "<<inputFile<<" for reading."<<std::endl;
        exit(1);
    }
    
    if ( (fp == NULL) && (strcmp(inputFile, "-") == 0) ) 
    {
        std::cerr<< PACKAGE_NAME<<" : [ERROR] Could not open stdin for reading."<<std::endl;
        exit(1);
    }
    return fp;
}

void RecursiveMkdir(std::string dir) 
{
    std::string tmp;
    size_t pos = 0;
    while ( std::string::npos != (pos = dir.find('/',pos+1)) ) {
        tmp = dir.substr(0,pos);
        mkdir(tmp.c_str(), S_IRWXU);
    }
    mkdir(dir.c_str(), S_IRWXU);
}


