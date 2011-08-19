
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
    
    std::string repeat, spacer, prev_spacer;
    repeat = spacer = prev_spacer = "";
    
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
        {   prev_spacer = spacer;
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
    int curr_repeat_start_index = mRepeats[i];
    int curr_repeat_end_index = curr_repeat_end_index + mRepeatLength;
    return mSequence.substr(curr_repeat_end_index, (curr_repeat_end_index - curr_repeat_start_index));
}


std::string Crispr::spacerStringAt(int i)
{
    int curr_repeat_end_index = mRepeats[i] + mRepeatLength - 1;
    int next_repeat_start_index = mRepeats[i + 1];
    int curr_spacer_start_index = curr_repeat_end_index + 1;
    int curr_spacer_end_index = next_repeat_start_index - 1;
    
    return mSequence.substr(curr_spacer_start_index, (curr_spacer_end_index - curr_spacer_start_index));
}

int Crispr::averageSpacerLength()
{
    int sum = 0;
    for (int i = 0; i < numSpacers(); i++)
    {
        sum += (int)spacerStringAt(i).length();
    }
    return sum/numSpacers();
}

int Crispr::averageRepeatLength()
{
    int sum = 0;
    for (int i = 0; i < numRepeats(); i++)
    {
        sum += (int)repeatStringAt(i).length();
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
}

int Crispr::getActualRepeatLength( int searchWindowLength, int minSpacerLength)
{
    int num_repeats = (int)(mRepeats.size() - 1);
    int first_repeat_start_index = mRepeats.front();
    int last_repeat_start_index = mRepeats.back();
    
    int shortest_repeat_spacing = mRepeats[1] - mRepeats[0];
    for (int i = 0; i < mRepeats.size() - 1; i++)
    {
        int curr_repeat_index = mRepeats[i];
        int next_repeat_index = mRepeats[i + 1];
        int curr_repeat_spacing = next_repeat_index - curr_repeat_index;
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }
    int sequence_length = (int)mSequence.length();
    
    int right_extension_length = searchWindowLength;
    int max_right_extension_length = shortest_repeat_spacing - minSpacerLength;
    
    
    int curr_repeat_start_index;
    std::string curr_repeat;
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    bool done = false;
    

    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
    while (!done && (right_extension_length <= max_right_extension_length) && (last_repeat_start_index + right_extension_length < sequence_length))
    {
        for (int k = 0; k < num_repeats; k++ )
        {
            curr_repeat_start_index = mRepeats.at(k);
            curr_repeat = mSequence.substr(curr_repeat_start_index, right_extension_length);
            char lastChar = curr_repeat.at(curr_repeat.length() - 1);
            
            if (lastChar == 'A')   char_count_A++;
            if (lastChar == 'C')   char_count_C++;
            if (lastChar == 'T')   char_count_T++;
            if (lastChar == 'G')   char_count_G++;
        }
        double percentA = (double)char_count_A/num_repeats;
        double percentC = (double)char_count_C/num_repeats;
        double percentT = (double)char_count_T/num_repeats;
        double percentG = (double)char_count_G/num_repeats;
        
        if ( (percentA >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentC >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentT >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentG >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            right_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            done = true;
        }
    }
    right_extension_length--;

    
    int left_extension_length = 0;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    done = false;
    
    int max_left_extension_length = shortest_repeat_spacing - minSpacerLength - right_extension_length;
    
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (!done && (left_extension_length <= max_left_extension_length) && (first_repeat_start_index - left_extension_length >= 0) )
    {
        for (int k = 0; k < num_repeats; k++ )
        {
            curr_repeat_start_index = mRepeats.at(k);
            char firstChar = mSequence.at(curr_repeat_start_index - left_extension_length);
            
            if (firstChar == 'A')    char_count_A++;
            if (firstChar == 'C')    char_count_C++;
            if (firstChar == 'T')    char_count_T++;
            if (firstChar == 'G')    char_count_G++;
        }
        
        double percentA = (double)char_count_A/num_repeats;
        double percentC = (double)char_count_C/num_repeats;
        double percentT = (double)char_count_T/num_repeats;
        double percentG = (double)char_count_G/num_repeats;
        
        if ( (percentA >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentC >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentT >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) || (percentG >= CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            left_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            done = true;
        }
    }
    left_extension_length--;

    repeatListIterator repeat_iter = mRepeats.begin();
    
    while (repeat_iter != mRepeats.end()) 
    {
        *repeat_iter = *repeat_iter - left_extension_length;
        ++repeat_iter;
    }
    mRepeatLength = right_extension_length + left_extension_length;
    return (right_extension_length + left_extension_length);
}


void Crispr::trim( int minRepeatLength)
{
    int num_repeats = this->numRepeats();
    
    std::string curr_repeat;
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    bool done = false;
    
    //trim from right
    while (!done && (this->repeatLength() > minRepeatLength) )
    {
        for (int k = 0; k < this->numRepeats(); k++ )
        {
            curr_repeat = this->repeatStringAt(k);
            char lastChar = curr_repeat[curr_repeat.length() - 1];
            
            if (lastChar == 'A')   char_count_A++;
            if (lastChar == 'C')   char_count_C++;
            if (lastChar == 'T')   char_count_T++;
            if (lastChar == 'G')   char_count_G++;
        }
        
        double percentA = (double)char_count_A/num_repeats;
        double percentC = (double)char_count_C/num_repeats;
        double percentT = (double)char_count_T/num_repeats;
        double percentG = (double)char_count_G/num_repeats;
        
        if ( (percentA < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentC < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentT < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentG < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            this->setRepeatLength(this->repeatLength() - 1);
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            done = true;
        }
    }
    
    
    
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    done = false;
    
    //trim from left
    while (!done && (this->repeatLength() > minRepeatLength) )
    {
        for (int k = 0; k < this->numRepeats(); k++ )
        {
            curr_repeat = this->repeatStringAt(k);
            char firstChar = curr_repeat.at(0);
            
            if (firstChar == 'A')   char_count_A++;
            if (firstChar == 'C')   char_count_C++;
            if (firstChar == 'T')   char_count_T++;
            if (firstChar == 'G')   char_count_G++;
    }
        
        double percentA = (double)char_count_A/num_repeats;
        double percentC = (double)char_count_C/num_repeats;
        double percentT = (double)char_count_T/num_repeats;
        double percentG = (double)char_count_G/num_repeats;
        
        if ( (percentA < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentC < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentT < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) && (percentG < CRASS_DEF_TRIM_EXTEND_CONFIDENCE) )
        {
            for (int m = 0; m < num_repeats; m++ )
            {
                int new_value = this->repeatAt(m) + 1;
                this->setRepeatAt(new_value, m);
            }
            this->setRepeatLength(this->repeatLength() - 1);
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
        {
            done = true;
        }
    }
}

bool Crispr::hasSimilarlySizedSpacers(void)
{
    int initial_spacer_length =(int)this->spacerStringAt(0).length();
    int repeat_length = this->repeatLength();
    
    for (int i = 0 ; i < this->numSpacers(); i++)
    {
        int curr_spacer_length = (int)this->spacerStringAt(i).length();
        
        //checks that each spacer is of similar size to other spacers
        if ( (curr_spacer_length - initial_spacer_length) >  CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF )
        {
            return false;
        }
        
        //checks that each spacer is of similar size to the repeats
        if ( (curr_spacer_length - repeat_length) > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF)
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
    std::string first_repeat = this->repeatStringAt(0);
    std::string first_spacer = this->spacerStringAt(0);
    float max_length;
    float similarity;
    float edit_distance;
    if (mRepeats.size() >= 3)
    {
        int i = 0;
        while ( (i < this->numSpacers() - 1) )
        {
            if (i == 4)  //only check first 5 spacers
                return true;
            
            std::string curr_spacer = this->spacerStringAt(i);
            std::string next_spacer = this->spacerStringAt(i + 1);
            std::string curr_repeat = this->repeatStringAt(i);
            
            max_length = std::max((int)curr_spacer.length(), (int)next_spacer.length());
            
            edit_distance =  LevenstheinDistance(curr_spacer, next_spacer);
            similarity = 1.0 - (edit_distance/max_length);
            
            //spacers should be different
            if ( similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
            {
                return false;
            }
            
            max_length = std::max((int)curr_repeat.length(), (int)next_spacer.length());
            edit_distance =  LevenstheinDistance(curr_repeat, next_spacer);
            similarity = 1.0 - (edit_distance/max_length);
            //repeats should also be different from spacers, otherwise may be tandem repeat
            if (similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
            {
                return false;
            }
            i++;
        }
        
        max_length = std::max((int)this->repeatStringAt(i).length(), (int)this->spacerStringAt(i).length());
        edit_distance =  LevenstheinDistance(this->repeatStringAt(i), this->spacerStringAt(i));
        similarity = 1.0 - (edit_distance/max_length);
        //checks last repeat/spacer
        if ( similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY )
        {
            return false;
        }
        return true;
    }
    
    //we check that the spacer is different from the repeat
    else if (mRepeats.size() == 2)
    {
        if (first_spacer == "")
        {
            return false;
        }
        else
        {
            max_length = std::max((int)first_spacer.length(), (int)first_repeat.length());
            edit_distance =  LevenstheinDistance(first_spacer, first_repeat);
            similarity = 1.0 - (edit_distance/max_length);
            if(similarity > CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY)
            {
                return false;
            }
        }
    }
    
    else
    {
        return false;
    }
    return true;
}

