/*
 *  Genome.cpp is part of the CRisprASSembler project
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


#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Genome.h"
#include "crass_defines.h"
#include "Levensthein.h"
#include <sstream>
std::string Crispr::toString()
{
    std::stringstream str;
    
    std::string repeat, spacer, prevSpacer;
    repeat = spacer = prevSpacer = "";
    
    str << "POSITION\tREPEAT\t\t\t\tSPACER"<<std::endl;
    
    str <<"--------\t";
    
    for (int y = 0; y <  this->repeatLength(); y++) str <<"-";
    str<<"\t";
    
    for (int z = 0; z <  this->averageSpacerLength(); z++) str << "-";
    str<<std::endl;
    
    
    //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
    for (int m = 0; m < this->numRepeats(); m++)
    {   //repeat = getRepeat(m);
        str << (this->repeatAt(m) + 1) << "\t\t" << repeatStringAt(m) << "\t";
        
        // print spacer
        // because there are no spacers after the last repeat, we stop early (m < crisprIndexVector.size() - 1)
        if (m < numSpacers())
        {   prevSpacer = spacer;
            spacer = spacerStringAt(m);
            str << spacer;
            
            str <<"\t[ " << this->repeatStringAt(m).length() << ", " << this->spacerStringAt(m).length() << " ]";
            //str +="--[" + DNASequence.getSimilarity(repeatStringAt(m), spacerStringAt(m)) + "]";
            //str +="--[" + DNASequence.getSimilarity(spacer, prevSpacer) + "]";
            str <<std::endl;
            
        }
    }
    
    
    str <<std::endl<<"--------\t";
    
    for (int x = 0; x < this->repeatLength(); x++)
    {
        str << "-";
    }
    str <<"\t";
    
    for (int z = 0; z <  this->averageSpacerLength(); z++)
    {
        str << "-";
    }
    str <<std::endl;
    
    
    return str.str();
    
}

std::string Crispr::repeatStringAt(int i)
{
    int currRepeatStartIndex = mRepeats[i];
    int currRepeatEndIndex = currRepeatStartIndex + mRepeatLength;
    return mSequence.substr(currRepeatStartIndex, (currRepeatEndIndex - currRepeatStartIndex));
}


std::string Crispr::spacerStringAt(int i)
{
    int currRepeatEndIndex = mRepeats[i] + mRepeatLength - 1;
    int nextRepeatStartIndex = mRepeats[i + 1];
    int currSpacerStartIndex = currRepeatEndIndex + 1;
    int currSpacerEndIndex = nextRepeatStartIndex - 1;
    
    return mSequence.substr(currSpacerStartIndex, (currSpacerEndIndex - currSpacerStartIndex));
}

int Crispr::averageSpacerLength()
{
    int sum = 0;
    for (int i = 0; i < numSpacers(); i++)
    {
        sum += spacerStringAt(i).length();
    }
    return sum/numSpacers();
}

int Crispr::averageRepeatLength()
{
    int sum = 0;
    for (int i = 0; i < numRepeats(); i++)
    {
        sum += repeatStringAt(i).length();
    }
    return sum/numRepeats();
}


// no remove element equivelent in c++
// going to have to hack this
void Crispr::removeRepeat(int val)
{
    // load what we've got into a tmp vector
    repeatList tmp_rep_list = mRepeats;
    mRepeats.clear();
    for (int i = 0; i < tmp_rep_list.size(); ++i) 
    {
        // skip the element
        if (i == val) 
        {
            continue;
        }
        // and add the rest
        else
        {
            mRepeats.push_back(tmp_rep_list.at(i));
        }
    }
    //mRepeats.insert(mRepeats.begin(),tmp_rep_list.begin(),tmp_rep_list.end());
}

int Crispr::getActualRepeatLength( int searchWindowLength, int minSpacerLength)
{
    int numRepeats = (int)(mRepeats.size() - 1);
    int firstRepeatStartIndex = mRepeats.front();
    int lastRepeatStartIndex = mRepeats.back();
    
    int shortestRepeatSpacing = mRepeats[1] - mRepeats[0];
    for (int i = 0; i < mRepeats.size() - 1; i++)
    {
        int currRepeatIndex = mRepeats[i];
        int nextRepeatIndex = mRepeats[i + 1];
        int currRepeatSpacing = nextRepeatIndex - currRepeatIndex;
        if (currRepeatSpacing < shortestRepeatSpacing)
            shortestRepeatSpacing = currRepeatSpacing;
    }
    int sequenceLength = (int)mSequence.length();
    
    int rightExtensionLength = searchWindowLength;
    int maxRightExtensionLength = shortestRepeatSpacing - minSpacerLength;
    
    
    int currRepeatStartIndex;
    std::string currRepeat;
    int charCountA, charCountC, charCountT, charCountG;
    charCountA = charCountC = charCountT = charCountG = 0;
    bool done = false;
    

    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
    while (!done && (rightExtensionLength <= maxRightExtensionLength) && (lastRepeatStartIndex + rightExtensionLength < sequenceLength))
    {
        for (int k = 0; k < numRepeats; k++ )
        {
            currRepeatStartIndex = mRepeats.at(k);
            currRepeat = mSequence.substr(currRepeatStartIndex, rightExtensionLength);
            std::cout<<"current repeat: "<<currRepeat<<std::endl;
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
        
        if ( (percentA >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentC >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentT >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentG >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            rightExtensionLength++;
            charCountA = charCountC = charCountT = charCountG = 0;
        }
        else
        {
            done = true;
        }
    }
    rightExtensionLength--;

    
    int leftExtensionLength = 0;
    charCountA = charCountC = charCountT = charCountG = 0;
    done = false;
    
    int maxLeftExtensionLength = shortestRepeatSpacing - minSpacerLength - rightExtensionLength;
    
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (!done && (leftExtensionLength <= maxLeftExtensionLength) && (firstRepeatStartIndex - leftExtensionLength >= 0) )
    {
        for (int k = 0; k < numRepeats; k++ )
        {
            currRepeatStartIndex = mRepeats.at(k);
            char firstChar = mSequence.at(currRepeatStartIndex - leftExtensionLength);
            
            if (firstChar == 'A')    charCountA++;
            if (firstChar == 'C')    charCountC++;
            if (firstChar == 'T')    charCountT++;
            if (firstChar == 'G')    charCountG++;
        }
        
        double percentA = (double)charCountA/numRepeats;
        double percentC = (double)charCountC/numRepeats;
        double percentT = (double)charCountT/numRepeats;
        double percentG = (double)charCountG/numRepeats;
        
        if ( (percentA >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentC >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentT >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentG >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            leftExtensionLength++;
            charCountA = charCountC = charCountT = charCountG = 0;
        }
        else
        {
            done = true;
        }
    }
    leftExtensionLength--;

    std::vector<int>::iterator newPositions = mRepeats.begin();
    
    while (newPositions != mRepeats.end()) 
    {
        *newPositions = *newPositions - leftExtensionLength;
        ++newPositions;
    }
    mRepeatLength = rightExtensionLength + leftExtensionLength;
    return (rightExtensionLength + leftExtensionLength);
}


void Crispr::trim( int minRepeatLength)
{
    int numRepeats = this->numRepeats();
    
    std::string currRepeat;
    int charCountA, charCountC, charCountT, charCountG;
    charCountA = charCountC = charCountT = charCountG = 0;
    bool done = false;
    
    //trim from right
    while (!done && (this->repeatLength() > minRepeatLength) )
    {
        for (int k = 0; k < this->numRepeats(); k++ )
        {
            currRepeat = this->repeatStringAt(k);
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
        
        if ( (percentA < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentC < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentT < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentG < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            this->setRepeatLength(this->repeatLength() - 1);
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
    while (!done && (this->repeatLength() > minRepeatLength) )
    {
        for (int k = 0; k < this->numRepeats(); k++ )
        {
            currRepeat = this->repeatStringAt(k);
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
        
        if ( (percentA < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentC < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentT < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentG < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            for (int m = 0; m < numRepeats; m++ )
            {
                int newValue = this->repeatAt(m) + 1;
                this->setRepeatAt(newValue, m);
            }
            this->setRepeatLength(this->repeatLength() - 1);
            charCountA = charCountC = charCountT = charCountG = 0;
        }
        else
        {
            done = true;
        }
    }
    
    //return candidateCRISPR;
}

bool Crispr::hasSimilarlySizedSpacers(void)
{
    int initialSpacerLength =(int)this->spacerStringAt(0).length();
    int repeatLength = this->repeatLength();
    
    for (int i = 0 ; i < this->numSpacers(); i++)
    {
        int currSpacerLength = (int)this->spacerStringAt(i).length();
        
        //checks that each spacer is of similar size to other spacers
        if ( (currSpacerLength - initialSpacerLength) >  CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF )
        {
            return false;
        }
        
        //checks that each spacer is of similar size to the repeats
        if ( (currSpacerLength - repeatLength) > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF)
        {
            return false;
        }
        
    }
    return true;
}


//checks first five spacers
bool Crispr::hasNonRepeatingSpacers(void)
{
    //assumes at least two elements
    std::string firstRepeat = this->repeatStringAt(0);
    std::string firstSpacer = this->spacerStringAt(0);
    float max_length;
    float similarity;
    float edit_distance;
    if ((mRepeats.size() - 1) >= 3)
    {
        int i = 0;
        while ( (i < this->numSpacers() - 1) )
        {
            if (i == 4)  //only check first 5 spacers
                return true;
            
            std::string currSpacer = this->spacerStringAt(i);
            std::string nextSpacer = this->spacerStringAt(i + 1);
            std::string currRepeat = this->repeatStringAt(i);
            
            max_length = std::max((int)currSpacer.length(), (int)nextSpacer.length());
            
            edit_distance =  Levensthein_distance(currSpacer, nextSpacer);
            similarity = 1.0 - (edit_distance/max_length);
            
            //spacers should be different
            if ( similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
                return false;
            
            max_length = std::max((int)currRepeat.length(), (int)nextSpacer.length());
            edit_distance =  Levensthein_distance(currRepeat, nextSpacer);
            similarity = 1.0 - (edit_distance/max_length);
            //repeats should also be different from spacers, otherwise may be tandem repeat
            if (similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
                return false;
            
            i++;
        }
        
        max_length = std::max((int)this->repeatStringAt(i).length(), (int)this->spacerStringAt(i).length());
        edit_distance =  Levensthein_distance(this->repeatStringAt(i), this->spacerStringAt(i));
        similarity = 1.0 - (edit_distance/max_length);
        //checks last repeat/spacer
        if ( similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
            return false;
        
        return true;
    }
    
    //we check that the spacer is different from the repeat
    else if ((mRepeats.size() - 1) == 2)
    {
        if (firstSpacer == "")
            return false;
        else
        {
            max_length = std::max((int)firstSpacer.length(), (int)firstRepeat.length());
            edit_distance =  Levensthein_distance(firstSpacer, firstRepeat);
            similarity = 1.0 - (edit_distance/max_length);
            if(similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY)
                return false;
        }
    }
    
    else
    {
        return false;
    }
    return true;
}

