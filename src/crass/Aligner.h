/*
 *  Aligner.h is part of the CRisprASSembler project
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

#ifndef crass_Aligner_h
#define crass_Aligner_h
#include <vector>
#include <bitset>
#include <map>

#include "ksw.h"
#include "StringCheck.h"
#include "Types.h"
#include "crassDefines.h"


typedef std::bitset<3> AlignerFlag_t;


class Aligner 
{
    // named accessors for the bitset flags
    enum Flag_t {
        reversed = 0,
        failed = 1,
        score_equal = 2
    };
    
    
    const unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    
    
    // ASCII table that converts characters into the multiplier
    // used for finding the correct index in the coverage array
    const unsigned char CHAR_TO_INDEX[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
public:
    //int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
    Aligner(int length, ReadMap *wh_reads, StringCheck *wh_st, int gapo=5, int gape=2, int minsc=5, int xtra=KSW_XSTART): 
        AL_length(length),
        AL_consensus(length,'N'), 
        AL_conservation(length, 0.0f), 
        AL_coverage(length*4, 0),
        AL_gapOpening(gapo), 
        AL_gapExtension(gape), 
        AL_minAlignmentScore(minsc), 
        AL_xtra(xtra) {
        
            // assign workhorse variables
            mReads = wh_reads;
            mStringCheck = wh_st;
        
        // set up default parameters for ksw alignment
        int sa = 1, sb = 3, i, j, k;
        AL_masterDRProfile = 0;
        AL_slaveDR_Profile = 0;

        if (AL_minAlignmentScore > 0xffff) AL_minAlignmentScore = 0xffff;
        if (AL_minAlignmentScore > 0) AL_xtra |= KSW_XSUBO | AL_minAlignmentScore;
        // initialize scoring matrix
        for (i = k = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j)
                AL_scoringMatrix[k++] = i == j? sa : -sb;
            AL_scoringMatrix[k++] = 0; // ambiguous base
        }
        for (j = 0; j < 5; ++j) AL_scoringMatrix[k++] = 0;
    }
    
    
    ~Aligner(){
        if (AL_masterDR != NULL) {
            delete AL_masterDR;
            AL_masterDR = NULL;
        } 
    }
    
    inline StringToken getMasterDrToken(){return AL_masterDRToken;}
    
    void setMasterDR(std::string& master);
    
    void alignSlave(StringToken& slaveDRToken);

    // add in all of the reads for this group to the coverage array
    void generateConsensus();
    
    inline std::map<StringToken, int>::iterator offsetBegin(){return AL_Offsets.begin();}
    
    inline std::map<StringToken, int>::iterator offsetEnd(){return AL_Offsets.end();}
    
    inline std::map<StringToken, int>::iterator offsetFind(StringToken& token){return AL_Offsets.find(token);}
    
    inline int offset(StringToken& drToken){return AL_Offsets[drToken];}

    inline int getDRZoneStart(){return AL_ZoneStart;}
    
    inline int getDRZoneEnd(){return AL_ZoneEnd;}
    
    inline void setDRZoneStart(int i){AL_ZoneStart = i;}
    
    inline void setDRZoneEnd(int i){AL_ZoneEnd = i;}
    
    inline int coverageAt(int i, char c){return AL_coverage.at(CHAR_TO_INDEX[c] * i); }
    
    inline char consensusAt(int i){return AL_consensus.at(i);}
    
    inline float conservationAt(int i){return AL_consensus.at(i);}
    
    inline int depthAt(int i){return AL_coverage[i] + AL_coverage[2*i] + AL_coverage[3*i] + AL_coverage[4*i];}

private:
    // private methods
    //
    
    // call ksw alignment to determine the offset for this slave against the master
    int getOffsetAgainstMaster(std::string& slaveDR,
                               AlignerFlag_t& flags);

    // transform any sequence into the right form for ksw
    void prepareSequenceForAlignment(std::string& sequence, uint8_t *transformedSequence);

    // transform a slave DR into the right form for ksw in both orientations
    void prepareSlaveForAlignment(std::string& slaveDR,
                                  uint8_t *slaveTransformedForward, 
                                  uint8_t *slaveTransformedReverse);

    // transform the master DR into the right form for ksw
    inline void prepareMasterForAlignment(std::string& masterDR) {
        AL_masterDRLength = masterDR.length();
        AL_masterDR = new uint8_t[AL_masterDRLength+1];
        prepareSequenceForAlignment(masterDR, AL_masterDR);
    };
    

    void placeReadsInCoverageArray(StringToken& currentDRToken);
    
    void extendSlaveDR(std::string& slaveDR, std::string& extendedSlaveDR);

    void calculateDRZone();
    
    
    
    
    //Members
    
    // length of the arrays
    int AL_length;
    
    // Vectors to hold the alignment data
    std::vector<char> AL_consensus;
    std::vector<float> AL_conservation;
    std::vector<int> AL_coverage;
    
    // Storage of all the offsets against the master
    std::map<StringToken, int> AL_Offsets;

    // smith-waterman default parameters
    int AL_gapOpening;
    int AL_gapExtension;
    int AL_minAlignmentScore;
    int AL_xtra;
    int8_t AL_scoringMatrix[25];
    
    // master DR
    uint8_t *AL_masterDR;
    int AL_masterDRLength;
    StringToken AL_masterDRToken;
    
    // "Glue" between WorkHorse
    ReadMap * mReads;
    StringCheck * mStringCheck;
    int AL_ZoneStart;
    int AL_ZoneEnd;

    
};



#endif
