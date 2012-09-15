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
#include "ksw.h"

class Aligner 
{
    static unsigned char seq_nt4_table[256] = {
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
public:
    //int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
    Aligner(int length, int gapo, int gape, int minsc, int xtra): 
        AL_consensus(length,'N'), 
        AL_conservation(length, 0.0f), 
        AL_coverage(length*4, 0),
        AL_gapOpening(gapo), 
        AL_gapExtension(gape), 
        AL_minAlignmentScore(minsc), 
        AL_xtra(xtra) {
        
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
    ~Aligner(){}

private:
    // private methods
    //

    // transform any sequence into the right form for ksw
    void prepareSequenceForAlignment(std::string& sequence, uint8_t *transformedSequence);

    // transform a slave DR into the right form for ksw in both orientations
    void prepareSlaveForAlignment(std::string& slaveDR,
                                  uint8_t *slaveTransformedForward, 
                                  uint8_t *slaveTransformedReverse);

    // transform the master DR into the right form for ksw
    inline void prepareMasterForAlignment(std::string& masterDR,
                                          uint8_t *masterTransformed) {
        return prepareSequenceForAlignment(masterDR, masterTransformed);
    };
    
    // call ksw alignment to determine the offset for this slave against the master
    int getOffsetAgainstMaster(uint8_t *slaveDrForward, 
                               uint8_t *slaveDrReverse, 
                               uint8_t *masterDR, 
                               bool& reversed);


    // Vectors to hold the alignment data
    std::vector<char> AL_consensus;
    std::vector<float> AL_conservation;
    std::vector<int> AL_coverage;

    // smith-waterman default parameters
    int AL_gapOpening;
    int AL_gapExtension;
    int AL_minAlignmentScore;
    int AL_xtra;
    int8_t AL_scoringMatrix[25];
    kswq_t *AL_masterDRProfile;
    kswq_t *AL_slaveDRProfile;

    
};



#endif
