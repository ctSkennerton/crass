/*
 *  SeqUtils.cpp is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011, 2012 Connor Skennerton & Michael Imelfort. All rights reserved. 
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

#include "SeqUtils.h"
#include <cstring>

char crass::comp_tab[128] = {
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
    std::string seq2 = crass::reverseComplement(seq1);
    if (seq1 < seq2)
    {
        return seq1;
    }
    return seq2;
}

char ** crass::cutIntoKmers(const char * s, int k, int& merCount)
{
    int str_len = static_cast<int>(strlen(s));
    int off = str_len - k;
    merCount = off + 1;

    
    // STOLED FROM SaSSY!!!!
    // First we cut kmers from the sequence then we use these to
    // determine overlaps, finally we make edges
    //
    // When we cut kmers from a read it is like this...
    //
    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    // ------------------------------
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    // XXXXXXXXXXXXXXXXXXXXXXXXXX
    //
    // so we break the job into three parts
    //
    // XXXX|XXXXXXXXXXXXXXXXXXXXXX|XXXX
    // ----|----------------------|----
    // XXXX|XXXXXXXXXXXXXXXXXXXXXX|
    //  XXX|XXXXXXXXXXXXXXXXXXXXXX|X
    //   XX|XXXXXXXXXXXXXXXXXXXXXX|XX
    //    X|XXXXXXXXXXXXXXXXXXXXXX|XXX
    //     |XXXXXXXXXXXXXXXXXXXXXX|XXXX
    //
    // the first and last part may be a little slow but the middle part can fly through...
    //
    
    // make a 2d array for the kmers!
	char ** kmers = NULL;
	int * kmer_offsets = NULL;
    kmers = new char*[merCount];
	
    for(int i = 0; i < merCount; i++)
    {
        kmers[i] = new char [k+1];
    }
    // use these offsets when we cut kmers, they are a component of the algorithm
    kmer_offsets = new int[merCount];
    for(int i = 0; i < merCount; i++)
    {
        kmer_offsets[i] = i * -1; // Starts at [0, -1, -2, -3, -4, ...]
    }
	
    int pos_counter = 0;
    
    // a slow-ish first part
    while(pos_counter < k)
    {
        for(int j = 0; j < merCount; j++)
        {
            if(pos_counter >= j)
            {
                kmers[j][kmer_offsets[j]] = s[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    
    // this is the fast part of the loop
    while(pos_counter < off)
    {
        for(int j = 0; j < merCount; j++)
        {
            if(kmer_offsets[j] >= 0 && kmer_offsets[j] < k)
            {
                kmers[j][kmer_offsets[j]] = s[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    
    // an even slower ending
    while(pos_counter < str_len)
    {
        for(int j = 0; j < merCount; j++)
        {
            if(kmer_offsets[j] < k)
            {
                kmers[j][kmer_offsets[j]] = s[pos_counter];
            }
            kmer_offsets[j]++;
        }
        pos_counter++;
    }
    delete [] kmer_offsets;
    return kmers;
}

void crass::cutIntoKmers(std::string& s, int k, std::vector<std::string>& mers)
{
    int off = static_cast<int>(s.length()) - k;
    for (int i = 0; i < off; ++i) {
        mers.push_back(s.substr(i,k));
    }
}




