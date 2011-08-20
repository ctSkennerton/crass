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
#include <zlib.h>  
#include <fstream>
#include <fcntl.h>
#include <stdlib.h>


// local includes
#include "libcrispr.h"
#include "LoggerSimp.h"
#include "crass_defines.h"
#include "WuManber.h"
#include "libbitap.h"
#include "bm.h"
#include "SeqUtils.h"
#include "Levensthein.h"
#include "Genome.h"

//**************************************
// DirectRepeat
//**************************************
DirectRepeat::DirectRepeat() {
    //-----
    // Constructor!
    //
    DR_Sequence = "";
    DR_MatchSequence = "";
    DR_Spacer = "";
    DR_Length = 0;
    DR_MatchStartPos = 0;
    DR_MatchEndPos = 0;
    DR_StartPos = 0;
    DR_EndPos = 0;
    DR_NumMismatches = 0;
}

void DirectRepeat::reset(void)
{
    //-----
    // DER
    //
    DR_Sequence = "";
    DR_MatchSequence = "";
    DR_Spacer = "";
    DR_Length = 0;
    DR_MatchStartPos = 0;
    DR_MatchEndPos = 0;
    DR_StartPos = 0;
    DR_EndPos = 0;
    DR_NumMismatches = 0;
}

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
        int seq_length = (int)strlen(seq->seq.s);//read.length();
        total_base += seq_length;
        read_counter++;
        if(read_counter > 100)
            break;
    }
    gzclose(fp);
    if((total_base / read_counter) > CRASS_DEF_READ_LENGTH_CUTOFF)
    {
        return LONG_READ;
    }
    return SHORT_READ;
}


// CRT search

float crtSearchFastqFile(const char *inputFastq, const options& opts, ReadMap * mReads, StringCheck * mStringCheck)
{
    logInfo("Long reads algorithm selected", 1);
    //-----
    // Code lifted from CRT, ported by connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    Crispr * candidate_crispr = new Crispr();
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        std::string read = seq->seq.s;
        std::string read_header = seq->name.s;
        int seq_length = (int)read.length() - 1;
        
        candidate_crispr->setSequence(read);
        
        total_base += seq_length;

        int actualRepeatLength;
        
        std::string pattern;
        
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
            pattern = read.substr(j, opts.searchWindowLength);
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
                actualRepeatLength = candidate_crispr->getActualRepeatLength(opts.searchWindowLength, opts.lowSpacerSize);
                if ( (actualRepeatLength >= opts.lowDRsize) && (actualRepeatLength <= opts.highDRsize) )
                {
                    if (candidate_crispr->hasNonRepeatingSpacers())
                    {
                        if (candidate_crispr->hasSimilarlySizedSpacers())
                        {
                            checkFlank(leftSide, candidate_crispr, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            checkFlank(rightSide, candidate_crispr, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            candidate_crispr->trim(opts.lowDRsize);

                            ReadHolder * tmp_holder = new ReadHolder();
                            Crispr::repeatListIterator rl_iter = candidate_crispr->mRepeats.begin();
                            while(rl_iter != candidate_crispr->mRepeats.end())
                            {
                                tmp_holder->RH_StartStops.push_back((*rl_iter));
                                //TODO -1 ?
                                tmp_holder->RH_StartStops.push_back((*rl_iter) + candidate_crispr->repeatLength());
                                rl_iter++;
                            }
                            //std::cout<<std::endl;
                            addReadHolder(mReads, mStringCheck, tmp_holder, read_header, read);
                            break;
                            //j = searchEnd + 1;
                        }
                        else
                        {
                            //candidate_crispr->superClear();
                            break;
                        }
                    }
                    else
                    {
                        //candidate_crispr->superClear();
                        break;
                    }
                 }
                else
                {
                    //candidate_crispr->superClear();
                    break;
                }
            }
        }
        candidate_crispr->superClear();
        read_counter++;
    }
    
    delete candidate_crispr;
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<inputFastq, 1);
    return total_base / read_counter;
}

// boyer moore functions
float bmSearchFastqFile(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    gzFile fp = getFileHandle(inputFastq);
    kseq_t *seq;
    int l, match_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    DirectRepeat dr_match;
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        std::string read = seq->seq.s;
        std::string read_header = seq->name.s;
        int seq_length = (int)read.length() - 1;
        int search_end = seq_length - opts.lowDRsize;
        
                
        logInfo("read counter: "<<match_counter, 8);
        
        total_base += seq_length;
        
        bool match_found = false;
        // create the read holder
        ReadHolder * tmp_holder = new ReadHolder();
        
        // boyer-moore search
        for (int start = 0; start < search_end; start++)
        {
            int search_begin = start + opts.lowDRsize + opts.lowSpacerSize;
            
            logInfo("search begin: " << search_begin << " search end: " << search_end, 8);
            
            if (search_begin >= search_end ) break;
            
            std::string query_word = read.substr(start, opts.lowDRsize);
            std::string subject_word = read.substr(search_begin);
            
            logInfo("query: " << query_word << " subject: " << subject_word, 10);
            
            int match_start_pos = PatternMatcher::bmpSearch( subject_word, query_word );
            logInfo("bm return: " << match_start_pos, 8);
            
            if (match_start_pos > -1) 
            {
                int end_pos = match_start_pos + search_begin + opts.lowDRsize;
                int matched_end_pos = start + opts.lowDRsize;
                
                dr_match.DR_StartPos = match_start_pos + search_begin;
                //dr_match.DR_StartList.push_back(MatchStartPos + search_begin);
                
                dr_match.DR_MatchStartPos = start;
                //dr_match.DR_StartList.push_back(start);
                
                dr_match.DR_MatchEndPos = matched_end_pos;
                dr_match.DR_EndPos = end_pos;
                
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                ++end_pos;
                ++matched_end_pos;
                if (end_pos <= seq_length) 
                {
                    // read through the subsuquent bases untill they don't match
                    logInfo("Read: "<<read_header<<" Len: "<<read.length(), 10);
                    while (read.at(matched_end_pos) == read.at(end_pos)) 
                    {
                        logInfo("Match end pos: "<<matched_end_pos<<" end pos: "<<end_pos, 10);
                        logInfo(read[matched_end_pos] << " : " << matched_end_pos << " == " << read[end_pos] << " : " << end_pos, 10);

                        
                        dr_match.DR_MatchEndPos = matched_end_pos;
                        dr_match.DR_EndPos = end_pos;
                        ++end_pos;
                        ++matched_end_pos;
                        if (end_pos > seq_length) break;
                    }
                    
                    if (cutDirectRepeatSequence(dr_match, opts, read))
                    {
                        patternsHash[dr_match.DR_MatchSequence] = true;

                        tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchStartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_MatchEndPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                        tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                        match_found = true;
                    } 
                    else 
                    {    
                        // increment the for loop so that it begins at the start of the 
                        // matched kmer/direct repeat
                        // minus 1 cause it will be incremented again at the top of the for loop
                        match_found = false;
                    }
                    start = dr_match.DR_StartPos - 1;
                    dr_match.reset();
                }
            }
        }
        if (match_found)
        {
            readsFound[read_header] = true;
            addReadHolder(mReads, mStringCheck, tmp_holder, read_header, read);
        }
        else
        {
            delete tmp_holder;
        }

        match_counter++;
    }
    
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<inputFastq, 1);
    return total_base / match_counter;
}


void findSingletons(const char *inputFastq, const options &opts, lookupTable &patternsHash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    logInfo("Beginning multipattern matcher: " << mReads->size(), 1);
    std::vector<std::string> patterns;
    
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
    
        DirectRepeat dr_match;        
        
        //initialize with an impossible number
        int start_pos = -1;
        
        std::string read = seq->seq.s;

        dr_match.DR_Sequence = search.Search(strlen(seq->seq.s), seq->seq.s, patterns, start_pos);
        dr_match.DR_StartPos = start_pos;
        
        if (start_pos != -1)
        {
            std::string header = seq->name.s;
            if (readsFound.find(header) == readsFound.end())
            //if (!(keyExists(readsFound, header)))
            {
                dr_match.DR_EndPos = start_pos + (int)dr_match.DR_Sequence.length();
                
                // create the read holder
                ReadHolder * tmp_holder = new ReadHolder;
                tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                addReadHolder(mReads, mStringCheck, tmp_holder, header, read);
            }
        }
    }
    logInfo("finished multi pattern matcher: " << mReads->size(), 1);
}

//**************************************
// kmer operators
//**************************************

bool cutDirectRepeatSequence(DirectRepeat& directRepeatMatch, const options& opts, std::string& read)
{
    directRepeatMatch.DR_Length = directRepeatMatch.DR_EndPos - directRepeatMatch.DR_StartPos;
    directRepeatMatch.DR_Sequence = read.substr(directRepeatMatch.DR_StartPos, (directRepeatMatch.DR_EndPos - directRepeatMatch.DR_StartPos + 1));
    directRepeatMatch.DR_MatchSequence = read.substr(directRepeatMatch.DR_MatchStartPos, (directRepeatMatch.DR_MatchEndPos - directRepeatMatch.DR_MatchStartPos + 1));
    directRepeatMatch.DR_Spacer = read.substr(directRepeatMatch.DR_MatchEndPos, (directRepeatMatch.DR_StartPos - directRepeatMatch.DR_MatchEndPos));

    if (!(checkDRAndSpacerLength(opts, directRepeatMatch)) || isLowComplexity(directRepeatMatch) || isSpacerAndDirectRepeatSimilar(directRepeatMatch))
    {
        return false;
    }

    // if the length of both spacer and direct repeat are okay cut the subsequences
    else
    {
        directRepeatMatch.DR_Sequence = read.substr(directRepeatMatch.DR_StartPos, (directRepeatMatch.DR_EndPos - directRepeatMatch.DR_StartPos + 1));
        directRepeatMatch.DR_MatchSequence = read.substr(directRepeatMatch.DR_MatchStartPos, (directRepeatMatch.DR_MatchEndPos - directRepeatMatch.DR_MatchStartPos + 1));
        directRepeatMatch.DR_Spacer = read.substr(directRepeatMatch.DR_MatchEndPos, (directRepeatMatch.DR_StartPos - directRepeatMatch.DR_MatchEndPos));
    } 
    return true;
}

bool checkDRAndSpacerLength(const options& opts, DirectRepeat& directRepeatMatch)
{
    logInfo("dr end pos: "<<directRepeatMatch.DR_EndPos<<" dr start pos: "<<directRepeatMatch.DR_StartPos, 10);
    
    int spacer_length = (int)directRepeatMatch.DR_Spacer.length();
    logInfo("DR len: "<<directRepeatMatch.DR_Length<<" SP length: "<<spacer_length, 10);
    // check if the direct repeat is in the right size range
    if ((directRepeatMatch.DR_Length < opts.lowDRsize) or (directRepeatMatch.DR_Length > opts.highDRsize)) 
    {
        return false;
    }
    
    // check if the spacer is in the right size range
    if ((spacer_length < opts.lowSpacerSize) or (spacer_length > opts.highSpacerSize))
    { 
        return false;
    }

    return true; 
}

bool isLowComplexity(DirectRepeat& directRepeatMatch)
{
    int cCount = 0;
    int gCount = 0;
    int aCount = 0;
    int tCount = 0;
    
    float aPercent;
    float cPercent;
    float gPercetn;
    float tPercent;

    std::string::iterator dr_iter = directRepeatMatch.DR_Sequence.begin();
    //int i = dr_match.DR_StartPos;
    while (dr_iter != directRepeatMatch.DR_Sequence.end()) 
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
            default: break;
        }
        dr_iter++;
    }
    aPercent = aCount/directRepeatMatch.DR_Length;
    tPercent = tCount/directRepeatMatch.DR_Length;
    gPercetn = gCount/directRepeatMatch.DR_Length;
    cPercent = cCount/directRepeatMatch.DR_Length;
    
    if (aPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || tPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || gPercetn > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || cPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        logInfo("Direct repeat has more than "<<CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD<<" of bases the same", 6);
        return true;   
    }
    return false;
}

bool isSpacerAndDirectRepeatSimilar(DirectRepeat& directRepeatMatch)
{
    float max_length = std::max((int)directRepeatMatch.DR_Spacer.length(), directRepeatMatch.DR_Length);
    
    float edit_distance = LevenstheinDistance(directRepeatMatch.DR_Sequence, directRepeatMatch.DR_Spacer);
    float similarity = 1.0 - (edit_distance/max_length);
    logInfo("similarity between spacer: "<<directRepeatMatch.DR_Spacer<<" and direct repeat: "<<directRepeatMatch.DR_Sequence<<" is: "<<similarity, 6);
    if (similarity > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        return true;
    }
    return false;
}

//**************************************
// transform read to DRlowlexi
//**************************************

std::string DRLowLexi(std::string& matchedRead, ReadHolder * tmpReadholder)
{
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    // make sure that tere is 4 elements in the array, if not you can only cut one
    if (tmpReadholder->RH_StartStops.size() == 2)
    {
        tmp_dr = matchedRead.substr(tmpReadholder->RH_StartStops.at(0), (tmpReadholder->RH_StartStops.at(1) - tmpReadholder->RH_StartStops.at(0) + 1));
        rev_comp = reverseComplement(tmp_dr);
    }
    else
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // if they both are then just take whichever is longer
        if (tmpReadholder->RH_StartStops.front() == 0 && tmpReadholder->RH_StartStops.back() == tmpReadholder->RH_Seq.length())
        {
            int lenA = tmpReadholder->RH_StartStops.at(1) - tmpReadholder->RH_StartStops.at(0);
            int lenB = tmpReadholder->RH_StartStops.at(3) - tmpReadholder->RH_StartStops.at(2);
            
            if (lenA > lenB)
            {
                tmp_dr = matchedRead.substr(tmpReadholder->RH_StartStops.at(0), (tmpReadholder->RH_StartStops.at(1) - tmpReadholder->RH_StartStops.at(0) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = matchedRead.substr(tmpReadholder->RH_StartStops.at(2), (tmpReadholder->RH_StartStops.at(3) - tmpReadholder->RH_StartStops.at(2) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
        }
        // take the second
        else if (tmpReadholder->RH_StartStops.front() == 0)
        {
            tmp_dr = matchedRead.substr(tmpReadholder->RH_StartStops.at(2), (tmpReadholder->RH_StartStops.at(3) - tmpReadholder->RH_StartStops.at(2) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
        // take the first
        else 
        {
            tmp_dr = matchedRead.substr(tmpReadholder->RH_StartStops.at(0), (tmpReadholder->RH_StartStops.at(1) - tmpReadholder->RH_StartStops.at(0) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
    }

    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        tmpReadholder->RH_WasLowLexi = true;
        tmpReadholder->RH_Seq = matchedRead;
        logInfo("DR in low lexi"<<endl<<tmpReadholder->RH_Seq, 9);
        return tmp_dr;
    }
    else
    {
        tmpReadholder->RH_Seq = reverseComplement(matchedRead);
        tmpReadholder->reverseStartStops();
        tmpReadholder->RH_WasLowLexi = false;
        logInfo("DR not in low lexi"<<endl<<tmpReadholder->RH_Seq, 9);
        return rev_comp;
    }
}

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmpReadholder, std::string& readHeader, std::string& read)
{
    logInfo("Add (header): \t" << tmpReadholder, 9);
    
    //add the header for the matched read
    tmpReadholder->RH_Header = readHeader;
    std::string dr_lowlexi = DRLowLexi(read, tmpReadholder);
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
    //out = new PrintStream(outputFileStream);
    
    //CTS//
    //outputFileStream<<"Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)\n");
    //outputFileStream<<"\n");
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
    //std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
    int num_repeats = candidateCRISPR->numRepeats();
    //std::cout<<"NR: "<<num_repeats<<std::endl;
    int pattern_length = (int)pattern.length();
    
    int last_repeat_index = candidateCRISPR->repeatAt(num_repeats-1);
    //std::cout<<"LR: "<<last_repeat_index<<std::endl;

    int second_last_repeat_index = candidateCRISPR->repeatAt(num_repeats-2);
    //std::cout<<"SLR: "<<second_last_repeat_index<<std::endl;

    int repeat_spacing = last_repeat_index - second_last_repeat_index;
    //std::cout<<"RS: "<<repeat_spacing<<std::endl;

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
        //        std::cout<<beginSearch<<" "<<mSequenceLength<<std::endl;
        //System.outputFileStream<<"beginSearch " + beginSearch + "  " + "endSearch" + endSearch);
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
        //        std::cout<<pattern<<"\t"<<text<<std::endl;
        position = PatternMatcher::bmpSearch(text, pattern);
        //        std::cout<<"bm pos: "<<position<<std::endl;
        
        
        if (position >= 0)
        {
            //std::cout<<"NP: "<<position<<std::endl;
            candidateCRISPR->addRepeat(begin_search + position);
            //std::cout<<"NEW: "<<begin_search + position<<std::endl;
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
