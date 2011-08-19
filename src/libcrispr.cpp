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

READ_TYPE decideWhichSearch(const char *input_fastq)
{
    //-----
    // Wrapper used for searching reads for DRs
    // depending on the length of the read. this funciton may use the boyer moore algorithm
    // of the CRT search algorithm
    //
    gzFile fp = getFileHandle(input_fastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;

    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        //std::string read = seq->seq.s;
        //std::string read_header = seq->name.s;
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

float crtSearchFastqFile(const char *input_fastq, const options &opts, ReadMap * mReads, StringCheck * mStringCheck)
{
    logInfo("Long reads algorithm selected", 1);
    //-----
    // Code lifted from CRT, ported by connor and hacked by Mike.
    // Should do well at finding crisprs in long reads
    //
    gzFile fp = getFileHandle(input_fastq);
    kseq_t *seq;
    int l, read_counter = 0;
    unsigned int total_base = 0;
    // initialize seq
    seq = kseq_init(fp);
    
    // read sequence  
    Crispr * candidateCRISPR = new Crispr();
    while ( (l = kseq_read(seq)) >= 0 ) 
    {
        std::string read = seq->seq.s;
        std::string read_header = seq->name.s;
        int seq_length = (int)read.length();
        
        candidateCRISPR->setSequence(read);
        
        total_base += seq_length;

        int actualRepeatLength;
        
        std::string pattern;
        
        //the mumber of bases that can be skipped while we still guarantee that the entire search
        //window will at some point in its iteration thru the sequence will not miss a any repeat
        int skips = opts.lowDRsize - (2 * opts.searchWindowLength - 1);
        if (skips < 1)
            skips = 1;

        int searchEnd = seq_length - opts.highDRsize - opts.highSpacerSize - opts.searchWindowLength;
        for (int j = 0; j <= searchEnd; j = j + skips)
        {
            
            int beginSearch = j + opts.lowDRsize + opts.lowSpacerSize;
            int endSearch = j + opts.highDRsize + opts.highSpacerSize + opts.searchWindowLength;
            
            if (endSearch > seq_length)
                endSearch = seq_length;
            
            if (endSearch < beginSearch) //should never occur
                endSearch = beginSearch;
            
            std::string text = read.substr(beginSearch, (endSearch - beginSearch));
            pattern = read.substr(j, opts.searchWindowLength);
            //if pattern is found, add it to candidate list and scan right for additional similarly spaced repeats
            int patternInTextIndex = PatternMatcher::bmpSearch(text, pattern);
            if (patternInTextIndex >= 0)
            {
                candidateCRISPR->addRepeat(j);
                candidateCRISPR->addRepeat(beginSearch + patternInTextIndex);
                scanRight(candidateCRISPR, pattern, opts.lowSpacerSize, 24, seq_length, read);
            }
            if ( (candidateCRISPR->numRepeats() >= opts.minNumRepeats) ) //make sure minNumRepeats is always at least 2
            {
                actualRepeatLength = candidateCRISPR->getActualRepeatLength(opts.searchWindowLength, opts.lowSpacerSize);
                if ( (actualRepeatLength >= opts.lowDRsize) && (actualRepeatLength <= opts.highDRsize) )
                {

                    if (candidateCRISPR->hasNonRepeatingSpacers())
                    {
                        if (candidateCRISPR->hasSimilarlySizedSpacers())
                        {
                            checkFlank(leftSide, candidateCRISPR, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            checkFlank(rightSide, candidateCRISPR, opts.lowSpacerSize, CRASS_DEF_SCAN_LENGTH, CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY, CRASS_DEF_SCAN_CONFIDENCE, seq_length, read);
                            candidateCRISPR->trim(opts.lowDRsize);
                            
                            ReadHolder * tmp_holder = new ReadHolder();
                            Crispr::repeatListIterator rl_iter = (candidateCRISPR->repeats()).begin();
                            while(rl_iter != (candidateCRISPR->repeats()).end())
                            {
                                tmp_holder->RH_StartStops.push_back(*rl_iter);
                                //TODO -1 ?
                                tmp_holder->RH_StartStops.push_back(*rl_iter + candidateCRISPR->repeatLength());
                                rl_iter++;
                            }
                            
                            addReadHolder(mReads, mStringCheck, tmp_holder, read_header, read);
                            break;
                            //j = searchEnd + 1;
                        }
                    }
                }
            }
        }
        candidateCRISPR->superClear();
        read_counter++;
    }
    
    delete candidateCRISPR;
    kseq_destroy(seq); // destroy seq  
    gzclose(fp);       // close the file handler  
    logInfo("finished processing file:"<<input_fastq, 1);
    return total_base / read_counter;
}

// boyer moore functions
float bmSearchFastqFile(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    gzFile fp = getFileHandle(input_fastq);
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
        int seq_length = read.length() - 1;
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
            
            int MatchStartPos = PatternMatcher::bmpSearch( subject_word, query_word );
            logInfo("bm return: " << MatchStartPos, 8);
            
            if (MatchStartPos > -1) 
            {
                int EndPos = MatchStartPos + search_begin + opts.lowDRsize;
                int MatchEndPos = start + opts.lowDRsize;
                
                dr_match.DR_StartPos = MatchStartPos + search_begin;
                //dr_match.DR_StartList.push_back(MatchStartPos + search_begin);
                
                dr_match.DR_MatchStartPos = start;
                //dr_match.DR_StartList.push_back(start);
                
                dr_match.DR_MatchEndPos = MatchEndPos;
                dr_match.DR_EndPos = EndPos;
                
                // make sure that the kmer match is not already at the end of the read before incrementing
                // increment so we are looking at the next base after the match
                ++EndPos;
                ++MatchEndPos;
                if (EndPos <= seq_length) 
                {
                    // read through the subsuquent bases untill they don't match
                    logInfo("Read: "<<seq->name.s<<" Len: "<<read.length(), 10);
                    while (read.at(MatchEndPos) == read.at(EndPos)) 
                    {
                        logInfo("Match end pos: "<<MatchEndPos<<" end pos: "<<EndPos, 10);
                        logInfo(read.at(MatchEndPos) << " : " << MatchEndPos << " == " << read.at(EndPos) << " : " << EndPos, 10);

                        
                        dr_match.DR_MatchEndPos = MatchEndPos;
                        dr_match.DR_EndPos = EndPos;
                        ++EndPos;
                        ++MatchEndPos;
                        if (EndPos > seq_length) break;
                    }
                    
                    if (cutDirectRepeatSequence(dr_match, opts, read))
                    {
                        patterns_hash[dr_match.DR_MatchSequence] = true;

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
    logInfo("finished processing file:"<<input_fastq, 1);
    return total_base / match_counter;
}


void findSingletons(const char *input_fastq, const options &opts, lookupTable &patterns_hash, lookupTable &readsFound, ReadMap * mReads, StringCheck * mStringCheck)
{
    logInfo("Beginning multipattern matcher: " << mReads->size(), 1);
    std::vector<std::string> patterns;
    
    map2Vector(patterns_hash, patterns);
    
    if (patterns.size() == 0)
    {
        logError("No patterns in vector for multimatch");
    }
    
    gzFile fp = getFileHandle(input_fastq);
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
                dr_match.DR_EndPos = start_pos + dr_match.DR_Sequence.length();
                
                // create the read holder
                ReadHolder * tmp_holder = new ReadHolder;
                tmp_holder->RH_StartStops.push_back(dr_match.DR_StartPos);
                tmp_holder->RH_StartStops.push_back(dr_match.DR_EndPos);
                addReadHolder(mReads, mStringCheck, tmp_holder, seq->name.s, read);
            }
        }
    }
    logInfo("finished multi pattern matcher: " << mReads->size(), 1);
}

//**************************************
// kmer operators
//**************************************

bool cutDirectRepeatSequence(DirectRepeat &dr_match, const options &opts, string &read)
{
    dr_match.DR_Length = dr_match.DR_EndPos - dr_match.DR_StartPos;
    dr_match.DR_Sequence = read.substr(dr_match.DR_StartPos, (dr_match.DR_EndPos - dr_match.DR_StartPos + 1));
    dr_match.DR_MatchSequence = read.substr(dr_match.DR_MatchStartPos, (dr_match.DR_MatchEndPos - dr_match.DR_MatchStartPos + 1));
    dr_match.DR_Spacer = read.substr(dr_match.DR_MatchEndPos, (dr_match.DR_StartPos - dr_match.DR_MatchEndPos));

    if (!(checkDRAndSpacerLength(opts, dr_match)) || isLowComplexity(dr_match) || isSpacerAndDirectRepeatSimilar(dr_match))
    {
        return false;
    }

    // if the length of both spacer and direct repeat are okay cut the subsequences
    else
    {
        dr_match.DR_Sequence = read.substr(dr_match.DR_StartPos, (dr_match.DR_EndPos - dr_match.DR_StartPos + 1));
        dr_match.DR_MatchSequence = read.substr(dr_match.DR_MatchStartPos, (dr_match.DR_MatchEndPos - dr_match.DR_MatchStartPos + 1));
        dr_match.DR_Spacer = read.substr(dr_match.DR_MatchEndPos, (dr_match.DR_StartPos - dr_match.DR_MatchEndPos));
    } 
    return true;
}

bool checkDRAndSpacerLength(const options &opts, DirectRepeat &dr_match)
{
    logInfo("dr end pos: "<<dr_match.DR_EndPos<<" dr start pos: "<<dr_match.DR_StartPos, 10);
    
    int spacer_length = (int)dr_match.DR_Spacer.length();
    logInfo("DR len: "<<dr_match.DR_Length<<" SP length: "<<spacer_length, 10);
    // check if the direct repeat is in the right size range
    if ((dr_match.DR_Length < opts.lowDRsize) or (dr_match.DR_Length > opts.highDRsize)) 
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

bool isLowComplexity(DirectRepeat &dr_match)
{
    int cCount = 0;
    int gCount = 0;
    int aCount = 0;
    int tCount = 0;
    
    float aPercent;
    float cPercent;
    float gPercetn;
    float tPercent;

    std::string::iterator dr_iter = dr_match.DR_Sequence.begin();
    //int i = dr_match.DR_StartPos;
    while (dr_iter != dr_match.DR_Sequence.end()) 
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
    aPercent = aCount/dr_match.DR_Length;
    tPercent = tCount/dr_match.DR_Length;
    gPercetn = gCount/dr_match.DR_Length;
    cPercent = cCount/dr_match.DR_Length;
    
    if (aPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || tPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || gPercetn > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD || cPercent > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        logInfo("Direct repeat has more than "<<CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD<<" of bases the same", 6);
        return true;   
    }
    return false;
}

bool isSpacerAndDirectRepeatSimilar(DirectRepeat &dr_match)
{
    float max_length = std::max((int)dr_match.DR_Spacer.length(), dr_match.DR_Length);
    
    float edit_distance = Levensthein_distance(dr_match.DR_Sequence, dr_match.DR_Spacer);
    float similarity = 1.0 - (edit_distance/max_length);
    logInfo("similarity between spacer: "<<dr_match.DR_Spacer<<" and direct repeat: "<<dr_match.DR_Sequence<<" is: "<<similarity, 6);
    if (similarity > CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD)
    {
        return true;
    }
    return false;
}

//**************************************
// transform read to DRlowlexi
//**************************************

std::string DRLowLexi(std::string matchedRead, ReadHolder * tmp_holder)
{
    //-----
    // Orientate a READ based on low lexi of the interalised DR
    //
    
    std::string tmp_dr;
    std::string rev_comp;
    
    // make sure that tere is 4 elements in the array, if not you can only cut one
    if (tmp_holder->RH_StartStops.size() == 2)
    {
        tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
        rev_comp = reverseComplement(tmp_dr);
    }
    else
    {
        // choose the dr that is not a partial ( no start at 0 or end at length)
        
        // if they both are then just take whichever is longer
        if (tmp_holder->RH_StartStops.front() == 0 && tmp_holder->RH_StartStops.back() == tmp_holder->RH_Seq.length())
        {
            int lenA = tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0);
            int lenB = tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2);
            
            if (lenA > lenB)
            {
                tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
            else
            {
                tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(2), (tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2) + 1));
                rev_comp = reverseComplement(tmp_dr);
            }
        }
        // take the second
        else if (tmp_holder->RH_StartStops.front() == 0)
        {
            tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(2), (tmp_holder->RH_StartStops.at(3) - tmp_holder->RH_StartStops.at(2) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
        // take the first
        else 
        {
            tmp_dr = matchedRead.substr(tmp_holder->RH_StartStops.at(0), (tmp_holder->RH_StartStops.at(1) - tmp_holder->RH_StartStops.at(0) + 1));
            rev_comp = reverseComplement(tmp_dr);
        }
    }

    
    if (tmp_dr < rev_comp)
    {
        // the direct repeat is in it lowest lexicographical form
        tmp_holder->RH_WasLowLexi = true;
        tmp_holder->RH_Seq = matchedRead;
        logInfo("DR in low lexi"<<endl<<tmp_holder->RH_Seq, 9);
        return tmp_dr;
    }
    else
    {
        tmp_holder->RH_Seq = reverseComplement(matchedRead);
        tmp_holder->reverseStartStops();
        tmp_holder->RH_WasLowLexi = false;
        logInfo("DR not in low lexi"<<endl<<tmp_holder->RH_Seq, 9);
        return rev_comp;
    }
}

void addReadHolder(ReadMap * mReads, StringCheck * mStringCheck, ReadHolder * tmp_holder, std::string read_header, std::string read)
{
    logInfo("Add (header): \t" << read_header, 9);
    
    //add the header for the matched read
    tmp_holder->RH_Header = read_header;
    std::string dr_lowlexi = DRLowLexi(read, tmp_holder);
    StringToken st = mStringCheck->getToken(dr_lowlexi);
    if(0 == st)
    {
        // new guy
        st = mStringCheck->addString(dr_lowlexi);
        (*mReads)[st] = new ReadList();
    }
    (*mReads)[st]->push_back(tmp_holder);
}


//**************************************
// STL extensions
//**************************************

// turn our map into a vector using just the keys
void map2Vector(lookupTable &patterns_hash, std::vector<std::string> &patterns)
{
    
    lookupTable::iterator iter = patterns_hash.begin();
    while (iter != patterns_hash.end()) 
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

void printGenomeCrispr(std::vector<Crispr*>& CRISPRVector, options& opts, std::string& mHeader, std::string& mSequence, bool repeatsFound)
{
    
    std::ofstream outputFileStream;
    
    std::string outputFileName = opts.output_fastq + "crass_genome.out";
    
    std::cout<<"Writing results in file " << outputFileName <<std::endl;            
    outputFileStream.open(outputFileName.c_str());
    //out = new PrintStream(outputFileStream);
    
    //CTS//
    //outputFileStream<<"Sequence '" + sequence.getName() + "' (" + sequence.length() + " bp)\n");
    //outputFileStream<<"\n");
    outputFileStream<<"ORGANISM: " << mHeader <<std::endl; 
    outputFileStream<<"BASES: " << mSequence.length() << std::endl<< std::endl;
    if (repeatsFound)
    {
        
        std::string repeat, spacer, prevSpacer;
        repeat = spacer = prevSpacer = "";
        
        //add 1 to each position, to offset programming languagues that begin at 0 rather than 1
        std::vector<Crispr*>::iterator crispr_iter = CRISPRVector.begin();
        int i = 1;
        while (crispr_iter != CRISPRVector.end()) 
        {
            outputFileStream<<"CRISPR"<<i<<" RANGE: "<< (*crispr_iter)->start() + 1 <<"-"<<(*crispr_iter)->end() + 1<<std::endl;           
            outputFileStream<<(*crispr_iter)->toString();
            
            outputFileStream<<"Total Repeats: " << (*crispr_iter)->numRepeats() <<" Average Repeat Length: " <<(*crispr_iter)->averageRepeatLength();
            outputFileStream<<" Average Spacer Length: " <<  (*crispr_iter)->averageSpacerLength()<<std::endl;
            outputFileStream<<std::endl<<std::endl;
            crispr_iter++; i++;
        }
        
    }
    else
    {
        outputFileStream<<"No CRISPR elements were found."<<std::endl;
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


void checkFlank(int side, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double spacerToSpacerMaxSimilarity, double confidence, int sequenceLength, std::string& read)
{
    bool moreToSearch = true;
    
    while (moreToSearch)
    {
        int result = scan(side, candidateCRISPR, minSpacerLength, scanRange, confidence, sequenceLength, read);
        if (result > 0)  //if another repeat found on flank
        {
            if (side == leftSide)
                candidateCRISPR->insertRepeatAt(result, 0);
            else if (side == rightSide)
                candidateCRISPR->addRepeat(result);
        }
        else
            moreToSearch = false;
    }
    
}

/*
 scan to the right and left of the first and last repeat to see if there is a region
 that is similar to the repeats.  necessary in case we missed a repeat because of
 inexact matches or a result of one of the filters
 */
int scan(int side, Crispr * candidateCRISPR, int minSpacerLength, int scanRange, double confidence, int mSequenceLength, std::string& mSequence)
{
    int repeatSpacing1, repeatSpacing2, avgRepeatSpacing;
    int firstRepeatIndex, lastRepeatIndex, candidateRepeatIndex;
    std::string repeatString, candidateRepeatString, newCandidateRepeatString;
    
    int repeatLength = candidateCRISPR->repeatLength();
    int numRepeats = candidateCRISPR->numRepeats();
    
    firstRepeatIndex = candidateCRISPR->repeatAt(0);
    lastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-1);
    
    if (side == leftSide)
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
    
    if (side == leftSide)
    {
        if (end > scanLeftMaxEnd)
            end = scanLeftMaxEnd;
    }
    
    if (side == rightSide)
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
    if ( (begin + repeatLength) > mSequenceLength)
    {
        return 0;
    }
    if ( (end + repeatLength) > mSequenceLength)
    {
        end = mSequenceLength - repeatLength;
    }
    
    if ( begin >= end)
        return 0;
    /******************** end range checks ********************/
    
    int array[end - begin + 1];
    
    int index = 0;
    for (int i = begin; i <= end; i++)
    {
        candidateRepeatString = mSequence.substr(i, repeatLength);
        array[index] = getHammingDistance(repeatString, candidateRepeatString);
        index++;
    }
    
    //min(array) returns the index of the smallest value in array  in this case, it refers to
    //the candidate string theat is closest to the repeatString.  uses hamming distance as levenshteinDistance is not useful for this particular task
    int newCandidateRepeatIndex = begin + min(array);
    newCandidateRepeatString = mSequence.substr(newCandidateRepeatIndex, repeatLength);
    
    bool match = patternMatches(repeatString, newCandidateRepeatString, confidence);
    
    if (match)
    {
        return newCandidateRepeatIndex;
    }
    else
    {
        return 0;
    }
    
}

void scanRight(Crispr * candidateCRISPR, std::string pattern, int minSpacerLength, int scanRange, int mSequenceLength, std::string& mSequence)
{
    int numRepeats = candidateCRISPR->numRepeats();
    int patternLength = (int)pattern.length();
    
    int lastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-1);
    
    int secondToLastRepeatIndex = candidateCRISPR->repeatAt(numRepeats-2);
    int repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;
    
    int candidateRepeatIndex, beginSearch, endSearch, position;
    
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
        {
            beginSearch = scanRightMinBegin;
        }
        //        std::cout<<beginSearch<<" "<<mSequenceLength<<std::endl;
        //System.outputFileStream<<"beginSearch " + beginSearch + "  " + "endSearch" + endSearch);
        if (beginSearch > mSequenceLength - 1)
        {
            return;
        }
        if (endSearch > mSequenceLength)
        {
            endSearch = mSequenceLength;
        }
        
        if ( beginSearch >= endSearch)
        {
            return;
        }
        /******************** end range checks ********************/
        
        std::string text = mSequence.substr(beginSearch, (endSearch - beginSearch));
        //        std::cout<<pattern<<"\t"<<text<<std::endl;
        position = PatternMatcher::bmpSearch(text, pattern);
        //        std::cout<<"bm pos: "<<position<<std::endl;
        
        if (position >= 0)
        {
            candidateCRISPR->addRepeat(beginSearch + position);
            secondToLastRepeatIndex = lastRepeatIndex;
            lastRepeatIndex = beginSearch + position;
            repeatSpacing = lastRepeatIndex - secondToLastRepeatIndex;
            if (repeatSpacing < (minSpacerLength + patternLength))
            {
                moreToSearch = false;
            }
        }
        else
        {
            moreToSearch = false;
        }
    }
    return;
}


bool patternMatches(std::string pattern1, std::string pattern2, double confidence)
{
    float max_length = std::max((int)pattern1.length(), (int)pattern2.length());
    float edit_distance =  Levensthein_distance(pattern1, pattern2);
    float similarity = 1.0 - (edit_distance/max_length);
    if (similarity >= confidence)
        return true;
    else
        return false;
}

int min (int * array)
{
    int min = array[0];
    int minIndex = 0;
    // get the number of elements in the array
    int length = (sizeof(*array)/sizeof(array[0]));
    
    for (int i = 0; i < length; i++)
    {
        if (array[i] < min)
        {
            min = array[i];
            minIndex = i;
        }
    }
    return minIndex;
}


int getHammingDistance(std::string seq1, std::string seq2)
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
