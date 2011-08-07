// File: crass_defines.h
// Original Author: Connor Skennerton on 7/05/11
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// The one stop shop for all your global definition needs!
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
// --------------------------------------------------------------------

#ifndef __CRASSDEFINES_H
    #define __CRASSDEFINES_H

#include <string>
#include <getopt.h>
// --------------------------------------------------------------------
 // PROGRAM NAME AND VERISON INFORMATION
// --------------------------------------------------------------------
#define CRASS_DEF_PRG_NAME                      "crass"
#define CRASS_DEF_LONG_NAME                     "CRisprASSembler"
#define CRASS_DEF_VERSION                       "0.0.1"
#define CRASS_DEF_MAJOR_VERSION                 0
#define CRASS_DEF_MINOR_VERSION                 1
#define CRASS_DEF_REVISION                      0
// --------------------------------------------------------------------
 // STRING LENGTH / MISMATCH / CLUSTER SIZE PARAMETERS
// --------------------------------------------------------------------
#define CRASS_DEF_READ_LENGTH_CUTOFF            (200.0)             // The length of read where the multipattern matcher would become insignificant. 
                                                                    //   ie where it is very likely that there is 2 full repeat-spacer-repeat units
#define CRASS_DEF_MAX_PATTERN_LENGTH            1024                // How long can a DR be?
#define CRASS_DEF_MAX_LOST_SOULS_MISMATCHES     2
#define CRASS_DEF_MAX_CLUSTER_SIZE_FOR_SW       30                  // Maximum number of cluster reps we will all vs all sw for
#define CRASS_DEF_MIN_SW_ALIGNMENT_RATIO        (0.85)              // SW alignments need to be this percentage of the original query to be considered real
#define CRASS_DEF_SW_SEARCH_EXT                 8
#define CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD     0.75
#define CRASS_DEF_KMER_SIZE                     7
#define CRASS_DEF_NUM_KMERS_4_MODE              5
// --------------------------------------------------------------------
// GENOME ALGORITHM DEFINES
// --------------------------------------------------------------------
#define  CRASS_DEF_SPACER_TO_SPACER_MAX_SIMILARITY 0.62
#define CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF 12
#define CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF 30
#define CRASS_DEF_MINREPEATIDENTITY 0.75
#define CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH 6
#define CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH 9
#define CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH 8
// --------------------------------------------------------------------
  // FILE IO
// --------------------------------------------------------------------
#define CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH     1024
#define CRASS_DEF_DEF_KMER_LOOKUP_EXT           "crass_kmers.txt"
#define CRASS_DEF_DEF_PATTERN_LOOKUP_EXT        "crass_direct_repeats.txt"
#define CRASS_DEF_DEF_SPACER_LOOKUP_EXT         "crass_spacers.txt"
// --------------------------------------------------------------------
 // USER OPTION STRUCTURE -- TODO: REMOVE THIS!
// --------------------------------------------------------------------
#define CRASS_DEF_MAX_DELIM_LENGTH 10                               // delimiter used in stats report
typedef struct {
    int count;
    bool detect;
    int logger_level;
    bool illumina;
    bool fourFiveFour;
    bool sanger;
    bool report_stats;
    int lowDRsize;
    int highDRsize;
    int lowSpacerSize;
    int highSpacerSize;
    int max_mismatches;
    std::string output_fastq;
    char delim[CRASS_DEF_MAX_DELIM_LENGTH];                         // delimiter used in stats report
    char * pat_file;
    int kmer_size;
} options;


static struct option long_options [] = {

    {"minDR", required_argument, NULL, 'd'},
    {"maxDR", required_argument, NULL, 'D'},
    {"minSpacer", required_argument, NULL, 's'},
    {"maxSpacer", required_argument, NULL, 'S'},
    {"logLevel", required_argument, NULL, 'l'},
    {"maxMismatches", required_argument, NULL, 'm'},
    {"version", no_argument, NULL, 'V'},
    {"kmerCount", required_argument, NULL, 'k'},
    {"outDir", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"reportStats", no_argument, NULL, 'r'},
    {"dumpReads", no_argument, NULL, 0},
    {"454", no_argument, NULL, 0},
    {"illumina", no_argument, NULL, 0},
    {"sanger", no_argument, NULL, 0},
    {NULL, no_argument, NULL, 0}
};
typedef struct {
    //std::string inputFileName;
    std::string outputFileDir;    // option 'o'
    int minNumRepeats;            // option 'n'
    int minRepeatLength;          // option 'd'
    int maxRepeatLength;          // option 'D'
    int minSpacerLength;          // option 's'
    int maxSpacerLength;          // option 'S'
    int searchWindowLength;       // option 'w'
    int logLevel;                 // option 'l'
} genOptions;

static struct option gen_long_options [] = {
    
    {"minDR", required_argument, NULL, 'd'},
    {"maxDR", required_argument, NULL, 'D'},
    {"minSpacer", required_argument, NULL, 's'},
    {"maxSpacer", required_argument, NULL, 'S'},
    {"logLevel", required_argument, NULL, 'l'},
    {"minNumRepeats", required_argument, NULL, 'n'},
    {"version", no_argument, NULL, 'V'},
    {"windowLength", required_argument, NULL, 'w'},
    {"outDir", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    //{"reportStats", no_argument, NULL, 'r'},
    //{"dumpReads", no_argument, NULL, 0},
    //{"454", no_argument, NULL, 0},
    //{"illumina", no_argument, NULL, 0},
    //{"sanger", no_argument, NULL, 0},
    //{"genome", no_argument, NULL, 0},
    {NULL, no_argument, NULL, 0}
};




// --------------------------------------------------------------------
#endif // __CRASSDEFINES_H
