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

void Aligner::prepareSequenceForAlignment(std::string& sequence, uint8_t *transformedSequence) {

    size_t seq_length = sequence.length();
    int i, j;
    for (i = 0; i < seq_length; ++i) 
        transformedSequence[i] = Aligner::seq_nt4_table[(int)sequence[i]];
    
    // null terminate the sequences
    transformedSequence[seq_length] = '\0';

}

void Aligner::prepareSlaveForAlignment(std::string& slaveDR)
