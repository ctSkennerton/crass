/*
 *  Aligner.cpp is part of the CRisprASSembler project
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
#include <iostream>
#include "Aligner.h"
#include "LoggerSimp.h"
#include "SeqUtils.h"

void Aligner::prepareSequenceForAlignment(std::string& sequence, uint8_t *transformedSequence) {

    size_t seq_length = sequence.length();
    size_t i;
    for (i = 0; i < seq_length; ++i) 
        transformedSequence[i] = seq_nt4_table[(int)sequence[i]];
    
    // null terminate the sequences
    transformedSequence[seq_length] = '\0';

}

void Aligner::prepareSlaveForAlignment(std::string& slaveDR, 
                                       uint8_t *slaveTransformedForward, 
                                       uint8_t *slaveTransformedReverse) {
    
    prepareSequenceForAlignment(slaveDR, slaveTransformedForward);
    std::string revcomp_slave_dr = reverseComplement(slaveDR);
    prepareSequenceForAlignment(revcomp_slave_dr, slaveTransformedReverse);
}

int Aligner::getOffsetAgainstMaster(std::string& slaveDR, AlignerFlag_t& flags) {

    int slave_dr_length = static_cast<int>(slaveDR.length());
    uint8_t slave_dr_forward[slave_dr_length+1];
    uint8_t slave_dr_reverse[slave_dr_length+1];
    
    prepareSlaveForAlignment(slaveDR, slave_dr_forward, slave_dr_reverse);

    // query profile 
    kswq_t *slave_forward_query_profile = 0;
    kswq_t *slave_reverse_query_profile = 0;
    
    
    
    // alignment of slave against master
    kswr_t forward_return = ksw_align(slave_dr_length, 
                                      slave_dr_forward, 
                                      AL_masterDRLength, 
                                      AL_masterDR, 
                                      5, 
                                      AL_scoringMatrix, 
                                      AL_gapOpening, 
                                      AL_gapExtension, 
                                      AL_xtra, 
                                      &slave_forward_query_profile);
    
    
    kswr_t reverse_return = ksw_align(slave_dr_length, 
                                      slave_dr_reverse, 
                                      AL_masterDRLength, 
                                      AL_masterDR, 
                                      5, 
                                      AL_scoringMatrix, 
                                      AL_gapOpening, 
                                      AL_gapExtension, 
                                      AL_xtra, 
                                      &slave_reverse_query_profile);
    
    
    // free the query profile
    free(slave_forward_query_profile); 
    free(slave_reverse_query_profile);
    
    // figure out which alignment was better
    if (reverse_return.score == forward_return.score) {
        flags[score_equal] = true;
    }
        
    if (reverse_return.score > forward_return.score && reverse_return.score >= AL_minAlignmentScore) {
        flags[reversed] = true;
        //std::cout<<"R: "<< masterDR<<"\t"<< reverse_return.tb<<"\t"<< reverse_return.te+1<<"\t"<< rev_slave<<"\t"<< reverse_return.qb<<"\t"<< reverse_return.qe+1<<"\t"<< reverse_return.score<<"\t"<< reverse_return.score2<<"\t"<< reverse_return.te2<<"\t"<<reverse_return.tb - reverse_return.qb<<std::endl;
        return reverse_return.tb - reverse_return.qb;
    } else if (forward_return.score >= AL_minAlignmentScore) {
        //std::cout<<"F: "<< masterDR<<"\t"<< forward_return.tb<<"\t"<< forward_return.te+1<<"\t"<< slaveDR<<"\t"<< forward_return.qb<<"\t"<< forward_return.qe+1<<"\t"<< forward_return.score<<"\t"<< forward_return.score2<<"\t"<< forward_return.te2<<"\t"<<forward_return.tb - forward_return.qb<<std::endl;
        return forward_return.tb - forward_return.qb;
    } else {
        logWarn("@Alignment Warning: Slave Score Failure",4);
        logWarn("Cannot place slave: "<<slaveDR<<" ("<<mStringCheck->getToken(slaveDR)<<") in array", 4);
        logWarn("Master: "<<AL_masterDR, 4);
        logWarn("Forward score: "<<forward_return.score, 4);
        logWarn("Reverse score: "<<reverse_return.score, 4);
        logWarn("******", 4);
        flags[failed] = true;
        return 0;
    }
    
    return 0;
}

void Aligner::generateCoverage() {
    
}

void Aligner::placeReadsInCoverageArray(StringToken& currentDrToken) {
    // we need to correct for the fact that we may not be using the 0th kmer
    int positional_offset = AL_Offsets[currentDrToken];//(kmer_positions_DR_master)[0] - (kmer_positions_DR_master)[positioning_kmer_index] + (kmer_positions_ARRAY)[positioning_kmer_index];
    ReadListIterator read_iter = mReads->at(currentDrToken)->begin();
    int current_dr_length = static_cast<int>(mStringCheck->getString(currentDrToken).length());
    
    while (read_iter != mReads->at(currentDrToken)->end()) 
    {
        // don't care about partials
        int dr_start_index = 0;
        int dr_end_index = 1;
        while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) != (current_dr_length - 1))
        {
            dr_start_index += 2;
            dr_end_index += 2;
        } 
        // go through every full length DR in the read and place in the array
        do
        {
            if(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == (current_dr_length - 1))
            {
                // we need to find the first kmer which matches the mode.
                int this_read_start_pos = positional_offset - (*read_iter)->startStopsAt(dr_start_index);
                for(int i = 0; i < (int)(*read_iter)->getSeqLength(); i++)
                {
                    int index = -1;
                    switch((*read_iter)->getSeqCharAt(i))
                    {
                        case 'A':
                            index = 1;
                            break;
                        case 'C':
                            index = 2;
                            break;
                        case 'G':
                            index = 3;
                            break;
                        default:
                            index = 4;
                            break;
                    }
                    int index_b = i+this_read_start_pos; 
                    if((index_b) >= AL_length)
                    {
                        logError("***FATAL*** MEMORY CORRUPTION: The consensus/coverage arrays are too short");
                    }
                    if((index_b) < 0)
                    {
                        logError("***FATAL*** MEMORY CORRUPTION: index = "<< index_b<<" less than array begining");
                    }
                    
                    AL_coverage[index * index_b]++;
                }
            }
            // go onto the next DR
            dr_start_index += 2;
            dr_end_index += 2;
            
            // check that this makes sense
            if(dr_start_index >= (int)((*read_iter)->numRepeats()*2)) {
                break;
            }
            
        } while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) == (current_dr_length - 1));
        read_iter++;
    }

}


void Aligner::extendSlaveDR(std::string &slaveDR, std::string &extendedSlaveDR){
 
    StringToken token = mStringCheck->getToken(slaveDR);
    
    // go into the reads and get the sequence of the DR plus a few bases on either side
    ReadListIterator read_iter = mReads->at(token)->begin();
    while (read_iter != mReads->at(token)->end()) 
    {
        // don't care about partials
        int dr_start_index = 0;
        int dr_end_index = 1;
        
        // Find the DR which is the right DR length.
        // compensates for partial repeats
        while(((*read_iter)->startStopsAt(dr_end_index) - (*read_iter)->startStopsAt(dr_start_index)) != ((int)(slaveDR.length()) - 1))
        {
            dr_start_index += 2;
            dr_end_index += 2;
        }
        // check that the DR does not lie too close to the end of the read so that we can extend
        if((*read_iter)->startStopsAt(dr_start_index) - 2 < 0 || (*read_iter)->startStopsAt(dr_end_index) + 2 > (*read_iter)->getSeqLength()) {
            // go to the next read
            read_iter++;
            continue;
        } else {
            // substring the read to get the new length
            extendedSlaveDR = (*read_iter)->getSeq().substr((*read_iter)->startStopsAt(dr_start_index) - 2, slaveDR.length() + 4);
        }
    }
}
















