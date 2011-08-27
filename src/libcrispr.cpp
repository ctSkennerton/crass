// File: libcrispr.cpp
// Original Author: Michael Imelfort 2011  :)
// --------------------------------------------------------------------
//
// OVERVIEW:
// 
// This file wraps all the crispr associated functions we'll need.
// The "crispr toolbox"
//
// --------------------------------------------------------------------
//  Copyright  2011 Michael Imelfort and Connor Skennerton
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//
// system includes
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <zlib.h>  
#include <fstream>
#include <fcntl.h>
#include <stdlib.h>


// local includes
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "WuManber.h"
#include "bm.h"
#include "SeqUtils.h"
#include "Levensthein.h"
#include "Crispr.h"


/* 
declare the type of file handler and the read() function
as described here:
http://lh3lh3.users.sourceforge.net/parsefastq.shtml

THIS JUST DEFINES A BUNCH OF **templated** structs

*/
KSEQ_INIT(gzFile, gzread);  

//**************************************
// search functions
//**************************************

READ_TYPE decideWhichSearch(const char *inputFastq)
{
    //-----
    // Wrapper used for searching reads for DRs
    // depending on the length of the read. this funciton may use the boyer moore algorithm
    // of the CRT search algorithm
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;

    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        int seq_length = (int)strlen(seq->seq.s);
        total_base += seq_length;
        read_counter++;
        if(read_counter > 100)
            break;
    }
    gzclose(fp);
    
    logInfo("Average read length (of the first 100 reads): "<<(total_base / read_counter), 1);
    if((total_base / read_counter) > CRASS_DEF_READ_LENGTH_CUTOFF)
    {
        return LONG_READ;
    }
    return SHORT_READ;
}


// CRT search

void longReadSearch(const char *inputFastq, const options& opts, ReadMap * mReads, StringCheck * mStringCheck)
{
    //-----
    // Code lifted from CRT, ported by connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    bool match_found = false;
    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    Crispr * candidate_crispr = new Crispr();
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s,seq->name.s);
        
        if (opts.fourFiveFour) 
        {
            tmp_holder->encode();
        }
        
        std::string read = tmp_holder->seq();
        
        int seq_length = (int)read.length() - 1;
        
        candidate_crispr->setSequence(read);
        
        total_base += seq_length;

        
        
        //the mumber of bases that can be skipped while we still guarantee that the entire search
        //window will at some point in its iteration thru the sequence will not miss a any repeat
        int skips = opts.lowDRsize - (2 * opts.searchWindowLength - 1);
        if (skips < 1)
        {
            skips = 1;
        }

        int searchEnd = seq_length - opts.highDRsize - opts.highSpacerSize - opts.searchWindowLength;
        for (int j = 0; j <= searchEnd; j = j + skips)
        {
            
            int beginSearch = j + opts.lowDRsize + opts.lowSpacerSize;
            int endSearch = j + opts.highDRsize + opts.highSpacerSize + opts.searchWindowLength;
            
            if (endSearch > seq_length)
            {
                endSearch = seq_length;
            }
            //should never occur
            if (endSearch < beginSearch)
            {
                endSearch = beginSearch;
            }
            
            std::string text = read.substr(beginSearch, (endSearch - beginSearch));
            std::string pattern = read.substr(j, opts.searchWindowLength);
            //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
            int pattern_in_text_index = PatternMatcher::bmpSearch(text, pattern);
            if (pattern_in_text_index >= 0)
            {
                candidate_crispr->addRepeat(j);
                candidate_crispr->addRepeat(beginSearch + pattern_in_text_index);

                scanRight(candidate_crispr, pattern, opts.lowSpacerSize, 24, seq_length, read);
            }

            if ( (candidate_crispr->numRepeats() >= opts.minNumRepeats) ) //make sure minNumRepeats is always at least 2
            {
                int actual_repeat_length = candidate_crispr->getActualRepeatLength(opts.searchWindowLength, opts.lowSpacerSize);
                if ( (actual_repeat_length >= opts.lowDRsize) && (actual_repeat_length <= opts.highDRsize) )
                {
                    if (candidate_crispr->hasNonRepeatingSpacers())
                    {
                        if (candidate_crispr->hasSimilarlySizedSpacers())
                        {
                            checkFlank(leftSide, candidate_crispr, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            checkFlank(rightSide, candidate_crispr, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            candidate_crispr->trim(opts.lowDRsize);
                            logInfo("potential CRISPR containing read found: "<<tmp_holder->header(), 3);
                            
                            match_found = true;
                            
                            Crispr::repeatListIterator rl_iter = candidate_crispr->mRepeats.begin();
                            while(rl_iter != candidate_crispr->mRepeats.end())
                            {
                                tmp_holder->add((*rl_iter), ((*rl_iter) + candidate_crispr->repeatLength()));
                                logInfo("direct repeat start index: "<<*rl_iter, 4);

                                rl_iter++;
                            }
                            logInfo(tmp_holder->seq(), 4);

                            addReadHolder(mReads, mStringCheck, tmp_holder, read);
                            break;
                        }
                        else
                        {
                            break;
                        }
                    }
                    else
                    {
                        break;
                    }
                 }
                else
                {
                    break;
                }
            }
        }
        if (!match_found) 
        {
            delete tmp_holder;
        }
        candidate_crispr->superClear();
        read_counter++;
    }
    
    delete candidate_crispr;
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<inputFastq, 1);
}

// boyer moore functions
void shortReadSearch(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    Crispr * candidate_crispr = new Crispr();
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        // create the read holder
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
        
        
        if (opts.fourFiveFour) 
        {
            tmp_holder->encode();
        } 
        std::string read = tmp_holder->seq();
        //std::string read = seq->seq.s;
        //std::string read_header = seq->name.s;
        
        candidate_crispr->setSequence(read);

        int seq_length = (int)read.length() - 1;
        int search_end = seq_length - opts.lowDRsize;
        
        total_base += seq_length;
        
        bool match_found = false;

        
        // boyer-moore search
        for (int start = 0; start < search_end; start++)
        {
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            
            if (search_begin >= search_end ) break;
            
            int match_start_pos = PatternMatcher::bmpSearch( read.substr(search_begin), read.substr(start, opts.lowDRsize) );
            
            if (match_start_pos > -1) 
            {
                int end_pos = match_start_pos + search_begin + opts.lowDRsize;
                int matched_end_pos = start + opts.lowDRsize;
                
                candidate_crispr->addRepeat(start);
                candidate_crispr->addRepeat(match_start_pos + search_begin);
                candidate_crispr->setRepeatLength(matched_end_pos - match_start_pos);

                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                ++end_pos;
                if (end_pos <= seq_length) 
                {
                    // read through the subsuquent bases untill they don't match
                    int extenstion_length = 0;
                    while (read.at(matched_end_pos) == read.at(end_pos)) 
                    {
                        extenstion_length++;
                        ++end_pos;
                        if (end_pos > seq_length) break;

                    }
                    candidate_crispr->setRepeatLength(candidate_crispr->repeatLength() + extenstion_length);
                }
                
                int repeat_length = candidate_crispr->repeatLength();
                int spacer_length = candidate_crispr->averageSpacerLength();
                if ( (repeat_length >= opts.lowDRsize) && (repeat_length <= opts.highDRsize) && (spacer_length >= opts.lowSpacerSize) && (spacer_length <= opts.highSpacerSize) )
                {
                    if (candidate_crispr->repeatAndSpacerIsDifferent()) 
                    {
                        if (!(candidate_crispr->isRepeatLowComplexity())) 
                        {
                            match_found = true;
                            logInfo("potential CRISPR containing read found: "<<tmp_holder->header(), 3);

                            patternsHash[candidate_crispr->repeatStringAt(0)] = true;
                            Crispr::repeatListIterator rep_iter = candidate_crispr->mRepeats.begin();
                            while(rep_iter != candidate_crispr->mRepeats.end())
                            {
                                tmp_holder->add((*rep_iter),((*rep_iter) + candidate_crispr->repeatLength()));
                                logInfo("direct repeat start index: "<<*rep_iter, 4);
                                rep_iter++;
                            }
                            logInfo(read, 4);

                        }
                    }
                }

                start = candidate_crispr->end() - 1;
                candidate_crispr->clear();
            }
        }
        if (match_found)
        {
            readsFound[tmp_holder->header()] = true;
            addReadHolder(mReads, mStringCheck, tmp_holder, read);
        }
        else
        {
            delete tmp_holder;
        }
        candidate_crispr->superClear();
        read_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<inputFastq, 1);
}


void findSingletons(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    logInfo("Beginning multipattern matcher. So far " << mReads->size()<<" reads have been found", 1);
    std::vector<std::string> patterns;
    int old_number = (int)mReads->size();
    map2Vector(patternsHash, patterns);
    
    if (patterns.size() == 0)
    {
        logError("No patterns in vector for multimatch");
    }
    
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    //ReadMatch match_info;
    
    seq = kseq_init(fp);

    WuManber search;
    search.Initialize(patterns);
    
    int l;
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        ReadHolder * tmp_holder = new ReadHolder(seq->seq.s, seq->name.s);
        
        if (opts.fourFiveFour) 
        {
            tmp_holder->encode();
        }
        std::string read = tmp_holder->seq();
        //initialize with an impossible number
        int start_pos = -1;
        std::string found_repeat = search.Search(read.length(), read.c_str(), patterns, start_pos);
        
        
        if (start_pos != -1)
        {
            if (readsFound.find(tmp_holder->header()) == readsFound.end())
            {
                logInfo("new read recruited: "<<tmp_holder->header(), 3);
                logInfo(tmp_holder->seq(), 4);
                tmp_holder->add(start_pos, (start_pos + (int)found_repeat.length()));
                addReadHolder(mReads, mStringCheck, tmp_holder, read);
            }
        }
        else
        {
            delete tmp_holder;
        }
    }
    logInfo("finished multi pattern matcher. An extra " << mReads->size() - old_number<<" reads were recruited", 1);
}


//**************************************
// transform read to DRlowlexi
//**************************************

std::string DRLowLexi( ReadHolder * tmpReadholder, std::string& read)
{
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    // make sure that tere is 4 elements in the array, if not you can only cut one
    if (tmpReadholder->drListSize() == 2)
    {
        tmp_dr = read.substr(tmpReadholder->at(0), (tmpReadholder->at(1) - tmpReadholder->at(0) + 1));
        rev_comp = reverseComplement(tmp_dr);
    }
    else
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // if they both are then just take whichever is longer
        if (tmpReadholder->front() == 0 && tmpReadholder->back() == read.length())
        {
            int lenA = tmpReadholder->at(1) - tmpReadholder->at(0);
            int lenB = tmpReadholder->at(3) - tmpReadholder->at(2);
            
            if (lenA > lenB)
            {
                tmp_dr = read.substr(tmpReadholder->at(0), (tmpReadholder->at(1) - tmpReadholder->at(0) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = read.substr(tmpReadholder->at(2), (tmpReadholder->at(3) - tmpReadholder->at(2) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
        }
        // take the second
        else if (tmpReadholder->at(0) == 0)
        {
            tmp_dr = read.substr(tmpReadholder->at(2), (tmpReadholder->at(3) - tmpReadholder->at(2) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
        // take the first
        else 
        {
            tmp_dr = read.substr(tmpReadholder->at(0), (tmpReadholder->at(1) - tmpReadholder->at(0) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
    }

    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        tmpReadholder->isLowLexi(true);
        //tmpReadholder->RH_Seq = matchedRead;
        logInfo("DR in low lexi"<<endl<<read, 9);
        return tmp_dr;
    }
    else
    {
        tmpReadholder->seq(reverseComplement(tmpReadholder->seq()));
        tmpReadholder->reverseStartStops();
        tmpReadholder->isLowLexi(false);
        logInfo("DR not in low lexi"<<endl<<tmpReadholder->seq(), 9);
        return rev_comp;
    }
}

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmpReadholder, std::string&  read)
{
    logInfo("Add (header): \t" << tmpReadholder->header(), 9);
    
    //add the header for the matched read
    //tmpReadholder->RH_Header = readHeader;
    std::string dr_lowlexi = DRLowLexi( tmpReadholder, read);
    StringToken st = mStringCheck->getToken(dr_lowlexi);
    if(0 == st)
    {
        // new guy
        st = mStringCheck->addString(dr_lowlexi);
        (*mReads)[st] = new ReadList();
    }
    (*mReads)[st]->push_back(tmpReadholder);
}

//**************************************
// STL extensions
//**************************************

// turn our map into a vector using just the keys
void map2Vector(lookupTable& patternsHash, std::vector<std::string>& patterns)
{
    
    lookupTable::iterator iter = patternsHash.begin();
    while (iter != patternsHash.end()) 
    {
        patterns.push_back(iter->first);
        iter++;
    }
}


//*************************************
// Old genome finder stuff
//*************************************
//bool goGenomeFinder(std::vector<std::string>& inputFiles)
//{
//    std::vector<std::string>::iterator inputFileName = inputFiles.begin();
//    while (inputFileName != inputFiles.end()) 
//    {
//        gzFile fp = getFileHandle(inputFileName->c_str());
//        std::cout<<std::endl<<std::endl<<"Reading file "<< *inputFileName<<std::endl;
//        
//        kseq_t *seq;
//        seq = kseq_init(fp);
//        int l;
//        // read sequence  
//        while ( (l = kseq_read(seq)) >= 0 ) 
//        {
//            mSequence = seq->seq.s;
//            mHeader = seq->name.s;
//            mSequenceLength = (int)mSequence.length();
//            std::cout<<mSequenceLength<<std::endl;
//            findRepeats();            
//        }
//        inputFileName++;
//    }
//    
//    logInfo("Done!",1);
//    return true;
//}
//

void printGenomeCrispr(std::vector<Crispr*>& CRISPRVector, options& opts, std::string& readHeader, std::string& read, bool repeatsFound)
{
    
    std::ofstream output_File_stream;
    
    std::string output_file_name = opts.output_fastq + "crass_genome.out";
    
    std::cout<<"Writing results in file " << output_file_name <<std::endl;            
    output_File_stream.open(output_file_name.c_str());

    output_File_stream<<"ORGANISM: " << readHeader <<std::endl; 
    output_File_stream<<"BASES: " << read.length() << std::endl<< std::endl;
    if (repeatsFound)
    {
        
        std::string repeat, spacer, prevSpacer;
        repeat = spacer = prevSpacer = "";
        
        //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
        std::vector<Crispr*>::iterator crispr_iter = CRISPRVector.begin();
        int i = 1;
        while (crispr_iter != CRISPRVector.end()) 
        {
            output_File_stream<<"CRISPR"<<i<<" RANGE: "<< (*crispr_iter)->start() + 1 <<"-"<<(*crispr_iter)->end() + 1<<std::endl;           
            output_File_stream<<(*crispr_iter)->toString();
            
            output_File_stream<<"Total Repeats: " << (*crispr_iter)->numRepeats() <<" Average Repeat Length: " <<(*crispr_iter)->averageRepeatLength();
            output_File_stream<<" Average Spacer Length: " <<  (*crispr_iter)->averageSpacerLength()<<std::endl;
            output_File_stream<<std::endl<<std::endl;
            crispr_iter++; i++;
        }
        
    }
    else
    {
        output_File_stream<<"No CRISPR elements were found."<<std::endl;
    }
    
    
    std::vector<Crispr*>::iterator crispr_iter = CRISPRVector.begin();
    while (crispr_iter != CRISPRVector.end()) 
    {
        if ((*crispr_iter) != NULL) 
        {
            delete (*crispr_iter);
            *crispr_iter = NULL;
        }
        crispr_iter++;
    }
    
}


void checkFlank(side sT, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence, int sequenceLength, std::string& read)
{
    bool more_to_search = true;
    
    while (more_to_search)
    {
        int result = scan(sT, candidateCRISPR, minSpacerLength, scanRange, confidence, sequenceLength, read);
        if (result > 0)  //if another repeat found on flank
        {
            if (sT == leftSide)
            {
                candidateCRISPR->insertRepeatAt(result, 0);
            }
            else if (sT == rightSide)
            {
                candidateCRISPR->addRepeat(result);
            }
        }
        else
        {
            more_to_search = false;
        }
    }
    
}

/*
 scan to the right and left of the first and last repeat to see if there is a region
 that is similar to the repeats.  necessary in case we missed a repeat because of
 inexact matches or a result of one of the filters
 */
int scan(side sT, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double confidence, int readLength, std::string& read)
{
    int repeat_spacing_1, repeat_spacing_2, avg_repeat_spacing;
    int first_repeat_index, last_repeat_index, candidate_repeat_index;
    std::string repeat_string, candidate_repeat_string, new_candidate_repeat_string;
    
    int repeat_length = candidateCRISPR->repeatLength();
    int num_repeats = candidateCRISPR->numRepeats();
    
    first_repeat_index = candidateCRISPR->repeatAt(0);
    last_repeat_index = candidateCRISPR->repeatAt(num_repeats-1);
    
    if (sT == leftSide)
    {
        repeat_string = candidateCRISPR->repeatStringAt(0);
        repeat_spacing_1 = candidateCRISPR->repeatSpacing(0, 1);
        if (num_repeats >= 3)
        {
            repeat_spacing_2 = candidateCRISPR->repeatSpacing(1, 2);
            avg_repeat_spacing = (repeat_spacing_1 + repeat_spacing_2)/2;
        }
        else
            avg_repeat_spacing = repeat_spacing_1;
        
        candidate_repeat_index = first_repeat_index - avg_repeat_spacing;
    }
    
    else //if (side.equals("right"))
    {
        repeat_string = candidateCRISPR->repeatStringAt(num_repeats-1);
        repeat_spacing_1 = candidateCRISPR->repeatSpacing(num_repeats-2, num_repeats-1);
        if (num_repeats >= 3)
        {
            repeat_spacing_2 = candidateCRISPR->repeatSpacing(num_repeats-3, num_repeats-2);
            avg_repeat_spacing = (repeat_spacing_1 + repeat_spacing_2)/2;
        }
        else
            avg_repeat_spacing = repeat_spacing_1;
        
        candidate_repeat_index = last_repeat_index + avg_repeat_spacing;
    }
    
    int begin = candidate_repeat_index - scanRange;
    int end   = candidate_repeat_index + scanRange;
    
    /******************** range checks ********************/
    //check that we do not search too far within an existing repeat when scanning right and left
    int scanLeftMaxEnd    = first_repeat_index - repeat_length - minSpacerLength;
    int scanRightMinBegin = last_repeat_index + repeat_length + minSpacerLength;
    
    if (sT == leftSide)
    {
        if (end > scanLeftMaxEnd)
            end = scanLeftMaxEnd;
    }
    
    if (sT == rightSide)
    {
        if (begin < scanRightMinBegin)
            begin = scanRightMinBegin;
    }
    
    //out of bounds check for scanning left
    if ( (begin) < 0)
    {
        return 0;
    }
    //out of bounds check for scanning right
    if ( (begin + repeat_length) > readLength)
    {
        return 0;
    }
    if ( (end + repeat_length) > readLength)
    {
        end = readLength - repeat_length;
    }
    
    if ( begin >= end)
        return 0;
    /******************** end range checks ********************/
    
    int array[end - begin + 1];
    
    int index = 0;
    for (int i = begin; i <= end; i++)
    {
        candidate_repeat_string = read.substr(i, repeat_length);
        array[index] = getHammingDistance(repeat_string, candidate_repeat_string);
        index++;
    }
    
    //min(array) returns the index of the smallest value in array  in this case, it refers to
    //the candidate string theat is closest to the repeatString.  uses hamming distance as levenshteinDistance is not useful for this particular task
    int new_candidate_repeat_index = begin + min(array);
    new_candidate_repeat_string = read.substr(new_candidate_repeat_index, repeat_length);
    
    bool match = patternMatches(repeat_string, new_candidate_repeat_string, confidence);
    
    if (match)
    {
        return new_candidate_repeat_index;
    }
    else
    {
        return 0;
    }
    
}

int scanRight(Crispr * candidateCRISPR, std::string& pattern, int minSpacerLength, int scanRange, int readLength, std::string& read)
{
    int num_repeats = candidateCRISPR->numRepeats();
    int pattern_length = (int)pattern.length();
    
    int last_repeat_index = candidateCRISPR->repeatAt(num_repeats-1);

    int second_last_repeat_index = candidateCRISPR->repeatAt(num_repeats-2);

    int repeat_spacing = last_repeat_index - second_last_repeat_index;

    int candidate_repeat_index, begin_search, end_search, position;
    
    bool more_to_search = true;
    while (more_to_search)
    {
        candidate_repeat_index = last_repeat_index + repeat_spacing;
        begin_search = candidate_repeat_index - scanRange;
        end_search = candidate_repeat_index + pattern_length + scanRange;
        
        /******************** range checks ********************/
        //check that we do not search too far within an existing repeat when scanning right
        int scanRightMinBegin = last_repeat_index + pattern_length + minSpacerLength;
        
        if (begin_search < scanRightMinBegin)
        {
            begin_search = scanRightMinBegin;
        }
        if (begin_search > readLength - 1)
        {
            return readLength - 1;
        }
        if (end_search > readLength)
        {
            end_search = readLength;
        }
        
        if ( begin_search >= end_search)
        {
            return end_search;
        }
        /******************** end range checks ********************/
        
        std::string text = read.substr(begin_search, (end_search - begin_search));
        position = PatternMatcher::bmpSearch(text, pattern);
        
        
        if (position >= 0)
        {
            candidateCRISPR->addRepeat(begin_search + position);
            second_last_repeat_index = last_repeat_index;
            last_repeat_index = begin_search + position;
            repeat_spacing = last_repeat_index - second_last_repeat_index;
            if (repeat_spacing < (minSpacerLength + pattern_length))
            {
                more_to_search = false;
            }
        }
        else
        {
            more_to_search = false;
        }
    }

    return begin_search + position;
}


bool patternMatches(std::string& pattern1, std::string& pattern2, double confidence)
{
    float max_length = std::max((int)pattern1.length(), (int)pattern2.length());
    float edit_distance =  LevenstheinDistance(pattern1, pattern2);
    float similarity = 1.0 - (edit_distance/max_length);
    if (similarity >= confidence)
        return true;
    else
        return false;
}

int min (int * array)
{
    int min = array[0];
    int min_index = 0;
    // get the number of elements in the array
    int length = (sizeof(*array)/sizeof(array[0]));
    
    for (int i = 0; i < length; i++)
    {
        if (array[i] < min)
        {
            min = array[i];
            min_index = i;
        }
    }
    return min_index;
}


int getHammingDistance(std::string& seq1, std::string& seq2)
{
    int length = (int)seq1.length();
    int hamming_distance = 0;
    
    if (seq1.length() != seq2.length())
    {
        length = (int)std::min(seq1.length(), seq2.length());
        hamming_distance = (int)seq1.length() - (int)seq2.length();
    }
    
    for (int i =0; i < length; i++)
    {
        if ( seq1.at(i) != seq2.at(i))
        {
            hamming_distance++;
        }
    }
    
    return hamming_distance;
}
