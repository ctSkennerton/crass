
/*
 *  Crispr.cpp is part of the CRisprASSembler project
 *
 *  There are already a number of options out there so instead of 
 *  reinventing the wheel, I've decided to port much of the code from 
 *  CRT with slight modifications.
 *
 *  Bland et al. (2007) "CRISPR Recognition Tool (CRT): a tool for automatic 
 *  detection of clustered regularly interspaced palindromic repeats" BMC 
 *  Bioinformatics 8:209.
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




#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "LoggerSimp.h"
#include "Crispr.h"
#include "crass_defines.h"
#include "Levensthein.h"


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
    int curr_repeat_end_index = curr_repeat_start_index + mRepeatLength;
    return mSequence.substr(curr_repeat_start_index, (curr_repeat_end_index - curr_repeat_start_index));
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
    for (int i = 0; i < (int)(tmp_rep_list.size()); ++i) 
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

int Crispr::extendPreRepeat(int searchWindowLength, int minSpacerLength)
{
    //-----
    // Extend a preliminary repeat - return the final repeat size
    //
    // the number of repeats
    int num_repeats = (int)(mRepeats.size());
    mRepeatLength = searchWindowLength;
    int cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * num_repeats);
        
    // the index in the read of the first DR kmer
    int first_repeat_start_index = mRepeats.front();
    
    // the index in the read of the last DR kmer
    int last_repeat_start_index = mRepeats.back();
    
    // the length between the first two DR kmers
    int shortest_repeat_spacing = mRepeats[1] - mRepeats[0];
    
    // loop througth all remaining members of mRepeats
    int end_index = (int)mRepeats.size();
    for (int i = 2; i < end_index; i++)
    {
        // get the repeat spacing of this pair of DR kmers
        int curr_repeat_spacing = mRepeats[i] - mRepeats[i - 1];

        // if it is shorter than what we already have, make it the shortest
        if (curr_repeat_spacing < shortest_repeat_spacing)
        {
            shortest_repeat_spacing = curr_repeat_spacing;
        }
    }

    // don't search too far  
    int max_right_extension_length = shortest_repeat_spacing - minSpacerLength;

    // Sometimes we shouldn't use the far right DR. (it may lie too close to the end)
    int DR_index_end = num_repeats;
    int dist_to_end = (int)mSequence.length() - last_repeat_start_index - 1;
    if(dist_to_end < max_right_extension_length)
    {
        DR_index_end--;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
    }
    std::string curr_repeat;
    int char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;


    //(from the right side) extend the length of the repeat to the right as long as the last base of all repeats are at least threshold
    while (max_right_extension_length > 0)
    {
        for (int k = 0; k < DR_index_end; k++ )
        {
            // look at the character just past the end of the last repeat
            switch(mSequence.at(mRepeats.at(k) + mRepeatLength))
            {
                case 'A':
                    char_count_A++;
                    break;
                case 'C':
                    char_count_C++;
                    break;
                case 'G':
                    char_count_G++;
                    break;
                case 'T':
                    char_count_T++;
                    break;
            }
        }
        
        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            mRepeatLength++;
            max_right_extension_length--;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
            break;
    }
    
    char_count_A = char_count_C = char_count_T = char_count_G = 0;

    // again, not too far
    int left_extension_length = 0;
    int max_left_extension_length = shortest_repeat_spacing - minSpacerLength - mRepeatLength;
    
    // and again, we may not wat to use the first DR
    int DR_index_start = 0;
    if(max_left_extension_length > first_repeat_start_index)
    {
        DR_index_start++;
        cut_off = (int)(CRASS_DEF_TRIM_EXTEND_CONFIDENCE * (num_repeats - 1));
    }
    
    //(from the left side) extends the length of the repeat to the left as long as the first base of all repeats is at least threshold
    while (left_extension_length < max_left_extension_length)
    {
        for (int k = DR_index_start; k < num_repeats; k++ )
        {
            switch(mSequence.at(mRepeats.at(k) - left_extension_length - 1))
            {
                case 'A':
                    char_count_A++;
                    break;
                case 'C':
                    char_count_C++;
                    break;
                case 'G':
                    char_count_G++;
                    break;
                case 'T':
                    char_count_T++;
                    break;
            }
        }

        if ( (char_count_A > cut_off) || (char_count_C > cut_off) || (char_count_G > cut_off) || (char_count_T > cut_off) )
        {
            mRepeatLength++;
            left_extension_length++;
            char_count_A = char_count_C = char_count_T = char_count_G = 0;
        }
        else
            break;
    }

    repeatListIterator repeat_iter = mRepeats.begin();
    
    while (repeat_iter != mRepeats.end()) 
    {
        *repeat_iter -= left_extension_length;
        if(*repeat_iter < 0)
            *repeat_iter = 0;
        ++repeat_iter;
    }

    return mRepeatLength;
}


void Crispr::trim( int minRepeatLength)
{
    float num_repeats = this->numRepeats();
    std::string curr_repeat;
    float char_count_A, char_count_C, char_count_T, char_count_G;
    char_count_A = char_count_C = char_count_T = char_count_G = 0;
    bool done = false;
    //trim from right
    while (!done && (this->repeatLength() > minRepeatLength) )
    {
        for (int k = 0; k < this->numRepeats(); k++ )
        {
            curr_repeat = this->repeatStringAt(k);
            char lastChar = curr_repeat[curr_repeat.length()-1];
            if (lastChar == 'A')   char_count_A++;
            if (lastChar == 'C')   char_count_C++;
            if (lastChar == 'T')   char_count_T++;
            if (lastChar == 'G')   char_count_G++;
        }
        float percentA = char_count_A/num_repeats;
        float percentC = char_count_C/num_repeats;
        float percentT = char_count_T/num_repeats;
        float percentG = char_count_G/num_repeats;
        
        
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
    
        
        float percentA = char_count_A/num_repeats;
        float percentC = char_count_C/num_repeats;
        float percentT = char_count_T/num_repeats;
        float percentG = char_count_G/num_repeats;
        
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
        int spacer_spacing = curr_spacer_length - initial_spacer_length;
        if ( spacer_spacing >  CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF )
        {
            logInfo("\tFailed test 4a. Spacer length's are too different to each other: "<<spacer_spacing<<" > " <<CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF, 8);
            return false;
        }
        int repeat_to_spacer_len_diff = curr_spacer_length - repeat_length;

        //checks that each spacer is of similar size to the repeats
        if ( repeat_to_spacer_len_diff > CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF)
        {
            logInfo("\tFailed test 4b. Spacer and Repeat lengths differ too much: "<<repeat_to_spacer_len_diff<<" > "<<CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF, 8);
            return false;
        }
        
    }
    return true;
}


//checks first five spacers
//assumes at least two elements
bool Crispr::hasNonRepeatingSpacers(void)
{
    try {
        if (mRepeats.size() < 2) 
        {
            logError("The vector holding the repeat indexes has a size less than 2!");
            throw "The vector holding the repeat indexes has a size less than 2!";
        }
    } catch (char * c) {
        std::cerr<<c<<std::endl;
    }
    if (mRepeats.size() >= 3)
    {
        int i = 0;
        while ( (i < this->numSpacers() - 1) )
        {
            //only check first 5 spacers
            if (i == 4) 
            {
                return true;
            }
            if (!this->areSpacersAtPosDifferent(i, i + 1))
            {
                return false;
            }
            if (!this->repeatAndSpacerIsDifferent(i)) 
            {
                return false;
            }
            i++;
        }
        
        //checks last repeat/spacer
        if (!this->repeatAndSpacerIsDifferent(i)) 
        {
            return false;
        }
        return true;
    }
    
    //we check that the spacer is different from the repeat
    else if (mRepeats.size() == 2)
    {
        if (this->spacerStringAt(0) == "")
        {
            return false;
        }
        else
        {
            return this->repeatAndSpacerIsDifferent();
        }
    }
    else
    {
        return false;
    }
    return true;
}
bool Crispr::areSpacersAtPosDifferent(int i, int j)
{
    std::string curr_spacer = this->spacerStringAt(i);
    std::string next_spacer = this->spacerStringAt(j);
    
    logInfo (curr_spacer << " : " << next_spacer, 1);
    float max_length = std::max(curr_spacer.length(), next_spacer.length());
    
    float edit_distance =  LevenstheinDistance(curr_spacer, next_spacer);
    float similarity = 1.0 - (edit_distance/max_length);
    if ( similarity > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY )
    {
        logInfo("\tFailed test 3a. Spacers are too similar: "<<similarity<< " > " << CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
        return false;
    }
    return true;
}

bool Crispr::repeatAndSpacerIsDifferent(int i)
{

    //assumes at least a single repeat and spacer elements
    std::string first_repeat = this->repeatStringAt(i);
    std::string first_spacer = this->spacerStringAt(i);

    
    float max_length = std::max(first_spacer.length(), first_repeat.length());
    float edit_distance =  LevenstheinDistance(first_spacer, first_repeat);
    float similarity = 1.0 - (edit_distance/max_length);
    if(similarity > CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY)
    {
        logInfo("\tFailed test 3b. Repeat and Spacer is too similar: "<<similarity<< " > " << CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY, 8);
        return false;
    }
    return true;
}

bool Crispr::isRepeatLowComplexity()
{
    float cCount = 0;
    float gCount = 0;
    float aCount = 0;
    float tCount = 0;
    float nCount = 0;
    
    float aPercent;
    float cPercent;
    float gPercetn;
    float tPercent;
    float nPercent;
    
    std::string curr_repeat = this->repeatStringAt(0);
    float curr_repeat_length = (float)curr_repeat.length();
    std::string::iterator dr_iter = curr_repeat.begin();
    //int i = dr_match.DR_StartPos;
    while (dr_iter != curr_repeat.end()) 
    {
        switch (*dr_iter) 
        {
            case 'c':
            case 'C':
                cCount++; break;
            case 't': 
            case 'T':
                tCount++; break;
            case 'a':
            case 'A':
                aCount++; break;
            case 'g':
            case 'G':
                gCount++; break;
            default: nCount++; break;
        }
        dr_iter++;
    }
    aPercent = aCount/curr_repeat_length;
    tPercent = tCount/curr_repeat_length;
    gPercetn = gCount/curr_repeat_length;
    cPercent = cCount/curr_repeat_length;
    nPercent = nCount/curr_repeat_length;
    

    if (aPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD) return true; 
    else if (tPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD) return true; 
    else if (gPercetn > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD) return true;
    else if (cPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD) return true;
    else if (nPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD) return true;   
    return false;
}

