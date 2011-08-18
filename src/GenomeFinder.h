/*
 *  GenomeFinder.h is part of the crass project
 *  
 *  Created by Connor Skennerton on 6/08/11.
 *  Copyright 2011 Connor Skennerton. All rights reserved. 
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
#ifndef crass_GenomeFinder_h
#define crass_GenomeFinder_h

#include <string>
#include <vector>

#include "Genome.h"
#include "crass_defines.h"
#include "kseq.h"


class GenomeFinder
{
public:
    enum side{right, left};

    GenomeFinder(genOptions * opts)
    {
        mOpts = opts;
    }
    ~GenomeFinder(void){}
    
    bool goGenomeFinder(std::vector<std::string>& inputFiles);
private:
    bool  findRepeats(void);
    void  trim(Crispr * candidateCRISPR, int minRepeatLength);
    void  checkFlank(int side, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence);
    int   scan(int side, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double confidence);
    void  scanRight(Crispr * candidateCRISPR, std::string pattern, int minSpacerLength, int scanRange);
    bool  patternMatches(std::string pattern1, std::string pattern2, double confidence);
    int   min (int * array);
    int   getHammingDistance(std::string seq1, std::string seq2);

    // members
    genOptions * mOpts;
//    std::string inputFileName;
//    std::string outputFileName;
//    
//    int screenDisplay;
//    int minNumRepeats;
//    int minRepeatLength;
//    int maxRepeatLength;
//    int minSpacerLength;
//    int maxSpacerLength;
//    int searchWindowLength;
    std::string mHeader;
    std::string mSequence;
    int mSequenceLength;
};

#endif
