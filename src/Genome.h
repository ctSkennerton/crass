/*
 *  Genome.h is part of the CRisprASSembler project
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


/*
    This code represents the algorithm for finding CRISPRs in genomes. 
    There are already a number of options out there so instead of 
    reinventing the wheel, I've decided to port much of the code from 
    CRT with slight modifications.
    
    Bland et al. (2007) "CRISPR Recognition Tool (CRT): a tool for automatic 
    detection of clustered regularly interspaced palindromic repeats" BMC 
    Bioinformatics 8:209.
 
 */

#ifndef crass_Genome_h
#define crass_Genome_h

#include <vector>
#include <string>

class GenomeCrispr
{
    typedef std::vector<int> repeatList;
    typedef std::vector<int>::iterator repeatListIterator;
public:
    // constriuctors
    GenomeCrispr()
    {
        repeatList mRepeats;
        mRepeatLength = 0;
        //mDnaSequence = _sequence;
    }
    
    GenomeCrispr( repeatList _positions, int _length)
    {
        mRepeats = _positions;
        mRepeatLength = _length;
        //mDnaSequence = _sequence;
    }
    
    //destructor
    ~GenomeCrispr()
    {
        //delete repeats;
    }
    
    // Getters and Setters
    inline repeatList repeats()
    {
        return mRepeats;
    }
    
    inline int repeatLength()
    {
        return mRepeatLength;
    }
    
    inline void setRepeats(repeatList _repeats)
    {
        mRepeats = _repeats;
    }
    
    inline void setRepeatLength(int length)
    {
        mRepeatLength = length;
    }
    
    inline int repeatSpacing(int pos1, int pos2)
    {
        return (repeatAt(pos2) - repeatAt(pos1));
    }
    
    inline void addRepeat(int val)
    {
        mRepeats.push_back(val);
    }
    
    inline void insertRepeatAt(int val, int pos)
    {
        repeatListIterator iter = mRepeats.begin()+pos;
        mRepeats.insert(iter, val);
    }
    
    inline void setRepeatAt(int val, int pos)
    {
        mRepeats[pos] = val;
    }
    
    inline int repeatAt(int i)
    {
        return mRepeats.at(i);
    }
    
    inline int start()
    {
        return mRepeats.front();
    }
    
    inline int end()
    {
        int lastRepeatBegin = mRepeats.back();
        return lastRepeatBegin + mRepeatLength - 1;
    }
    
    inline int firstRepeat()
    {
        return mRepeats.front();
    }
    
    inline int lastRepeat()
    {
        return mRepeats.back();
    }
    
    inline int numRepeats()
    {
        return (int)mRepeats.size() - 1;
    }
    
    inline int numSpacers()
    {
        return numRepeats() - 1;
    }
    
    
    std::string repeatStringAt(int i);
    std::string spacerStringAt(int i);
    int averageSpacerLength(void);
    int averageRepeatLength(void);
    std::string toString(void);
    void removeRepeat(int val);
    bool hasSimilarlySizedSpacers(void);
    bool hasNonRepeatingSpacers(void);
    int getActualRepeatLength( int searchWindowLength, int minSpacerLength);
    void trim( int minRepeatLength);

    
    friend class GenomeFinder;
private:
    repeatList mRepeats;
    int mRepeatLength;
};
    
    



#endif
