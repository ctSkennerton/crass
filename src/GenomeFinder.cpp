/*
 *  GenomeFinder.cpp is part of the crass project
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

#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>
#include <fstream>

#include "Genome.h"
#include "GenomeFinder.h"
#include "LoggerSimp.h"
#include "kseq.h"
#include "SeqUtils.h"
#include "bm.h"
#include "crass_defines.h"
/* 
 declare the type of file handler and the read() function
 as described here:
 http://lh3lh3.users.sourceforge.net/parsefastq.shtml
 
 THIS JUST DEFINES A BUNCH OF **templated** structs
 
 */
KSEQ_INIT(gzFile, gzread);

bool GenomeFinder::goGenomeFinder(std::vector<std::string>& inputFiles)
{
    std::vector<std::string>::iterator inputFileName = inputFiles.begin();
    while (inputFileName != inputFiles.end()) 
    {
        gzFile fp = getFileHandle(inputFileName->c_str());
        std::cout<<std::endl<<std::endl<<"Reading file "<< *inputFileName<<std::endl;
        
        kseq_t *seq;
        seq = kseq_init(fp);
        int l;
        // read sequence  
        while ( (l = kseq_read(seq)) >= 0 ) 
        {
            mSequence = seq->seq.s;
            mHeader = seq->name.s;
            mSequenceLength = (int)mSequence.length();

            //CTS//
            //System.out.println("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)")
            try
            {
                findRepeats();
            }
            catch (int e)
            {
                std::cerr<<"Error processing input file "<< *inputFileName << ". Please, check contents."<<std::endl;
            }
            
        }

    }
        
    logInfo("Done!",1);
    return true;
}
    
    
bool GenomeFinder::findRepeats(void)
{
    
    std::vector<GenomeCrispr*> CRISPRVector;
    int actualRepeatLength;
    bool repeatsFound = false;
    
    GenomeCrispr * candidateCRISPR  = new GenomeCrispr();
    std::string pattern;
    

    
    //the mumber of bases that can be skipped while we still guarantee that the entire search
    //window will at some point in its iteration thru the sequence will not miss a any repeat
    int skips = mOpts->minRepeatLength - (2 * mOpts->searchWindowLength - 1);
    if (skips < 1)
        skips = 1;
    
    //CTS//
    //System.out.println("Searching for repeats...");
    //long repeatSearchStart = System.currentTimeMillis();
    
    //SearchUtil searchUtil = new SearchUtil();
    
    int searchEnd = mSequenceLength - mOpts->maxRepeatLength - mOpts->maxSpacerLength - mOpts->searchWindowLength;
    for (int j = 0; j <= searchEnd; j = j + skips)
    {
        //candidateCRISPR = new CRISPR();
        
        int beginSearch = j + mOpts->minRepeatLength + mOpts->minSpacerLength;
        int endSearch = j + mOpts->maxRepeatLength + mOpts->maxSpacerLength + mOpts->searchWindowLength;
        
        if (endSearch > mSequenceLength)
            endSearch = mSequenceLength;
        
        if (endSearch < beginSearch)  //should never occur
            endSearch = beginSearch;
        
        std::string text = mSequence.substr(beginSearch, endSearch);
        pattern = mSequence.substr(j, j + mOpts->searchWindowLength);
        
        //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
        int patternInTextIndex = PatternMatcher::bmpSearch(text, pattern);
        if (patternInTextIndex >= 0)
        {
            candidateCRISPR->addRepeat(j);
            candidateCRISPR->addRepeat(beginSearch + patternInTextIndex);
            this->scanRight(candidateCRISPR, pattern, mOpts->minSpacerLength, 24);
        }
        
        if ( (candidateCRISPR->numRepeats() >= mOpts->minNumRepeats) )  //make sure minNumRepeats is always at least 2
        {
            actualRepeatLength = candidateCRISPR->getActualRepeatLength(mOpts->searchWindowLength, mOpts->minSpacerLength);
            
            if ( (actualRepeatLength >= mOpts->minRepeatLength) && (actualRepeatLength <= mOpts->maxRepeatLength) )
            {   
                if (candidateCRISPR->hasNonRepeatingSpacers())
                {
                    if (candidateCRISPR->hasSimilarlySizedSpacers())
                    {
                         this->checkFlank("left", candidateCRISPR, mOpts->minSpacerLength, 30, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, .70);
                         this->checkFlank("right", candidateCRISPR, mOpts->minSpacerLength, 30, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, .70);
                         this->trim(candidateCRISPR, mOpts->minRepeatLength);
                         CRISPRVector.push_back(candidateCRISPR);
                         repeatsFound = true;
                        
                        //we may skip current CRISPR (assuming CRISPRs are not interleaved)
                        j = candidateCRISPR->end() + 1;
                    }
                }
            }
        }
    }
    
    //CTS//
    //long repeatSearchEnd = System.currentTimeMillis();
    //System.out.println("Time to search for repeats:  " + (repeatSearchEnd - repeatSearchStart) + " ms");
    //System.out.println(CRISPRVector.size() + " possible CRISPR(s) found" + "\n");
    
    
    // ********************** Display CRISPR elements ********************** //
    // ********************************************************************* //
    try
    {
        
        std::ofstream outputFileStream;

        std::string outputFileName = mOpts->outputFileDir + "a.out";
            
        std::cout<<"Writing results in file " << outputFileName <<std::endl;            
        outputFileStream.open(outputFileName.c_str());
            //out = new PrintStream(outputFileStream);
        
        //CTS//
        //out.print("Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)\n");
        //out.print("\n");
        
        if (repeatsFound)
        {
            outputFileStream<<">" << mHeader << "_" << mSequence.length() << "_";
            
            //out.print("\n");
            //int repeatLength, numRepeats, numSpacers;
            //CRISPR currCRISPR;
            
            std::string repeat, spacer, prevSpacer;
            repeat = spacer = prevSpacer = "";
            
            //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
            std::vector<GenomeCrispr*>::iterator crispr_iter = CRISPRVector.begin();
            while (crispr_iter != CRISPRVector.end()) 
            {
                //currCRISPR = (CRISPR)CRISPRVector.elementAt(k);
                outputFileStream<<"crispr_" << ((*crispr_iter)->start() + 1) << "-" <<  ((*crispr_iter)->end() + 1) << "_";
                outputFileStream<<"RN_" + (*crispr_iter)->numRepeats() << "_RL_" << (*crispr_iter)->averageRepeatLength() + "_SL_" <<  (*crispr_iter)->averageSpacerLength() <<std::endl;
                outputFileStream<<mSequence<<std::endl;
                outputFileStream<<"Total Repeats: " << (*crispr_iter)->numRepeats() <<" Repeat Length: " <<(*crispr_iter)->averageRepeatLength();
                outputFileStream<<" Average Spacer Length: " <<  (*crispr_iter)->averageSpacerLength()<<std::endl;
                outputFileStream<<(*crispr_iter)->toString();
                outputFileStream<<std::endl<<std::endl;
            }
            for (int k = 0; k < CRISPRVector.size(); k++)
            {

            }
            //CTS//
            //out.print("Time to find repeats: " + (repeatSearchEnd - repeatSearchStart) + " ms\n\n");
        }
        else
        {
            outputFileStream<<"No CRISPR elements were found."<<std::endl;
        }
    }
    catch (int e)   {   std::cerr<<"--Error writing to file-- "<<std::endl;   }
    
    delete candidateCRISPR;
    //delete [] CRISPRVector;
    
    return true;
}


void GenomeFinder::trim(GenomeCrispr * candidateCRISPR, int minRepeatLength)
    {
        int numRepeats = candidateCRISPR->numRepeats();
//        int left = (int)candidateCRISPR->start();
//        int right = (int)candidateCRISPR->end();
//        
        std::string currRepeat;
        int charCountA, charCountC, charCountT, charCountG;
        charCountA = charCountC = charCountT = charCountG = 0;
        bool done = false;
        
        //trim from right
        while (!done && (candidateCRISPR->repeatLength() > minRepeatLength) )
        {
            for (int k = 0; k < candidateCRISPR->numRepeats(); k++ )
            {
                currRepeat = candidateCRISPR->repeatStringAt(k);
                char lastChar = currRepeat.at(currRepeat.length() - 1);
                
                if (lastChar == 'A')   charCountA++;
                if (lastChar == 'C')   charCountC++;
                if (lastChar == 'T')   charCountT++;
                if (lastChar == 'G')   charCountG++;
            }
            
            double percentA = (double)charCountA/numRepeats;
            double percentC = (double)charCountC/numRepeats;
            double percentT = (double)charCountT/numRepeats;
            double percentG = (double)charCountG/numRepeats;
            
            if ( (percentA < .75) && (percentC < .75) && (percentT < .75) && (percentG < .75) )
            {
                candidateCRISPR->setRepeatLength(candidateCRISPR->repeatLength() - 1);
                charCountA = charCountC = charCountT = charCountG = 0;
            }
            else
            {
                done = true;
            }
        }
        
        
        
        charCountA = charCountC = charCountT = charCountG = 0;
        done = false;
        
        //trim from left
        while (!done && (candidateCRISPR->repeatLength() > minRepeatLength) )
        {
            for (int k = 0; k < candidateCRISPR->numRepeats(); k++ )
            {
                currRepeat = candidateCRISPR->repeatStringAt(k);
                char firstChar = currRepeat.at(0);
                
                if (firstChar == 'A')   charCountA++;
                if (firstChar == 'C')   charCountC++;
                if (firstChar == 'T')   charCountT++;
                if (firstChar == 'G')   charCountG++;
            }
            
            double percentA = (double)charCountA/numRepeats;
            double percentC = (double)charCountC/numRepeats;
            double percentT = (double)charCountT/numRepeats;
            double percentG = (double)charCountG/numRepeats;
            
            if ( (percentA < .75) && (percentC < .75) && (percentT < .75) && (percentG < .75) )
            {
                for (int m = 0; m < numRepeats; m++ )
                {
                    int newValue = candidateCRISPR->repeatAt(m) + 1;
                    candidateCRISPR->setRepeatAt(newValue, m);
                }
                candidateCRISPR->setRepeatLength(candidateCRISPR->repeatLength() - 1);
                charCountA = charCountC = charCountT = charCountG = 0;
            }
            else
            {
                done = true;
            }
        }
        
        //return candidateCRISPR;
    }
    
 
    
void GenomeFinder::checkFlank(std::string side, GenomeCrispr * candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence)
    {
        bool moreToSearch = true;
        
        while (moreToSearch)
        {
            int result = scan(side, candidateCRISPR, minSpacerLength, scanRange, confidence);
            if (result > 0)  //if another repeat found on flank
            {
                if (side == "left")
                    candidateCRISPR->insertRepeatAt(result, 0);
                else if (side == "right")
                    candidateCRISPR->addRepeat(result);
            }
            else
                moreToSearch = false;
        }
        
        //return candidateCRISPR;
    }
    
    /*
     scan to the right and left of the first and last repeat to see if there is a region
     that is similar to the repeats.  necessary in case we missed a repeat because of
     inexact matches or a result of one of the filters
     */
int GenomeFinder::scan(std::string side, GenomeCrispr * candidateCRISPR, int minSpacerLength, int scanRange, double confidence)
    {
        int repeatSpacing1, repeatSpacing2, avgRepeatSpacing;
        int firstRepeatIndex, lastRepeatIndex, candidateRepeatIndex;
        std::string repeatString, candidateRepeatString, newCandidateRepeatString;
        
        int repeatLength = candidateCRISPR->repeatLength();
        int numRepeats = candidateCRISPR->numRepeats();
        
        firstRepeatIndex = candidateCRISPR->repeatAt(0);
        lastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-1);
        
        if (side == "left")
        {
            repeatString = candidateCRISPR->repeatStringAt(0);
            repeatSpacing1 = candidateCRISPR->repeatSpacing(0, 1);
            if (numRepeats >= 3)
            {
                repeatSpacing2 = candidateCRISPR->repeatSpacing(1, 2);
                avgRepeatSpacing = (repeatSpacing1 + repeatSpacing2)/2;
            }
            else
                avgRepeatSpacing = repeatSpacing1;
            
            candidateRepeatIndex = firstRepeatIndex - avgRepeatSpacing;
        }
        
        else //if (side.equals("right"))
        {
            repeatString = candidateCRISPR->repeatStringAt(numRepeats-1);
            repeatSpacing1 = candidateCRISPR->repeatSpacing(numRepeats-2, numRepeats-1);
            if (numRepeats >= 3)
            {
                repeatSpacing2 = candidateCRISPR->repeatSpacing(numRepeats-3, numRepeats-2);
                avgRepeatSpacing = (repeatSpacing1 + repeatSpacing2)/2;
            }
            else
                avgRepeatSpacing = repeatSpacing1;
            
            candidateRepeatIndex = lastRepeatIndex + avgRepeatSpacing;
        }
        
        int begin = candidateRepeatIndex - scanRange;
        int end   = candidateRepeatIndex + scanRange;
        
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right and left
        int scanLeftMaxEnd    = firstRepeatIndex - repeatLength - minSpacerLength;
        int scanRightMinBegin = lastRepeatIndex + repeatLength + minSpacerLength;
        
        if (side == "left")
        {
            if (end > scanLeftMaxEnd)
                end = scanLeftMaxEnd;
        }
        
        if (side == "right")
        {
            if (begin < scanRightMinBegin)
                begin = scanRightMinBegin;
        }
        
        //out of bounds check for scanning left
        if ( (begin) < 0)
            return 0;
        
        //out of bounds check for scanning right
        if ( (begin + repeatLength) > mSequenceLength)
            return 0;
        if ( (end + repeatLength) > mSequenceLength)
            end = mSequenceLength - repeatLength;
        
        if ( begin >= end)
            return 0;
        /******************** end range checks ********************/
        
        int array[end - begin + 1];
        
        int index = 0;
        for (int i = begin; i <= end; i++)
        {
            candidateRepeatString = mSequence.substr(i, i + repeatLength);
            array[index] = this->getHammingDistance(repeatString, candidateRepeatString);
            index++;
        }
        
        //min(array) returns the index of the smallest value in array  in this case, it refers to
        //the candidate string theat is closest to the repeatString.  uses hamming distance as levenshteinDistance is not useful for this particular task
        int newCandidateRepeatIndex = begin + min(array);
        newCandidateRepeatString = mSequence.substr(newCandidateRepeatIndex, newCandidateRepeatIndex + repeatLength);
        
        bool match = patternMatches(repeatString, newCandidateRepeatString, confidence);
        
        if (match)
            return newCandidateRepeatIndex;
        else
            return 0;
        
    }
    
void GenomeFinder::scanRight(GenomeCrispr * candidateCRISPR, std::string pattern, int minSpacerLength, int scanRange)
{
    int numRepeats = candidateCRISPR->numRepeats();
    int patternLength = (int)pattern.length();
    
    int lastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-1);
    int secondToLastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-2);
    int repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;
    
    int candidateRepeatIndex, beginSearch, endSearch, position;
    std::string text = "";
    
    bool moreToSearch = true;
    while (moreToSearch)
    {
        candidateRepeatIndex = lastRepeatIndex + repeatSpacing;
        beginSearch = candidateRepeatIndex - scanRange;
        endSearch = candidateRepeatIndex + patternLength + scanRange;
        
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right
        int scanRightMinBegin = lastRepeatIndex + patternLength + minSpacerLength;
        
        if (beginSearch < scanRightMinBegin)
            beginSearch = scanRightMinBegin;
        
        //System.out.print("beginSearch " + beginSearch + "  " + "endSearch" + endSearch);
        if (beginSearch > mSequenceLength - 1)
            return;
        if (endSearch > mSequenceLength)
            endSearch = mSequenceLength;
        
        if ( beginSearch >= endSearch)
            return;
        /******************** end range checks ********************/
        
        text = mSequence.substr(beginSearch, endSearch);
        position = PatternMatcher::bmpSearch(text, pattern);
        
        if (position >= 0)
        {
            candidateCRISPR->addRepeat(beginSearch + position);
            secondToLastRepeatIndex = lastRepeatIndex;
            lastRepeatIndex = beginSearch + position;
            repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;
            if (repeatSpacing < (minSpacerLength + patternLength))
                moreToSearch = false;
        }
        else
            moreToSearch = false;
    }
    
    return;
}
    
    
bool GenomeFinder::patternMatches(std::string pattern1, std::string pattern2, double confidence)
    {
        double patternSimilarity = DNASequence.getSimilarity(pattern1, pattern2);
        if (patternSimilarity >= confidence)
            return true;
        else
            return false;
    }
    
int GenomeFinder::min (int * array)
{
    int min = array[0];
    int minIndex = 0;
    
    for (int i = 0; i < array.length; i++)
    {
        if (array[i] < min)
        {
            min = array[i];
            minIndex = i;
        }
    }
    return minIndex;
}


int GenomeFinder::getHammingDistance(std::string seq1, std::string seq2)
{
    int length = (int)seq1.length();
    int hammingDistance = 0;
    
    if (seq1.length() != seq2.length())
    {
        length = (int)std::min(seq1.length(), seq2.length());
        hammingDistance = (int)seq1.length() - (int)seq2.length();
    }
    
    for (int i =0; i < length; i++)
    {
        if ( seq1.at(i) != seq2.at(i))
        {
            hammingDistance++;
        }
    }
    
    return hammingDistance;
}

