
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
#include "Rainbow.h"
#include "config.h"
// --------------------------------------------------------------------
 // General Macros
// --------------------------------------------------------------------
#define isNotDecimal(number) (number > 1.0 || number < 0.0)
// --------------------------------------------------------------------
 // STRING LENGTH / MISMATCH / CLUSTER SIZE PARAMETERS
// --------------------------------------------------------------------
#define CRASS_DEF_MAX_CLUSTER_SIZE_FOR_SW       30                  // Maximum number of cluster reps we will all vs all sw for
#define CRASS_DEF_MIN_SW_ALIGNMENT_RATIO        (0.85)              // SW alignments need to be this percentage of the original query to be considered real
#define CRASS_DEF_SW_SEARCH_EXT                 8
#define CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD     (0.75)
#define CRASS_DEF_KMER_SIZE                     7
#define CRASS_DEF_K_CLUST_MIN                   12
#define CRASS_DEF_READ_COUNTER_LOGGER           100000
#define CRASS_DEF_MAX_READS_FOR_DECISION        1000
// --------------------------------------------------------------------
// HARD CODED PARAMS FOR FINDING TRUE DRs
// --------------------------------------------------------------------
#define CRASS_DEF_PERCENT_IN_ZONE_CUT_OFF       (0.85)              // amount that a DR must agrre with the existsing DR within a zone to be added
#define CRASS_DEF_NUM_KMERS_4_MODE              (5)                 // find the top XX occuring in the DR
#define CRASS_DEF_NUM_KMERS_4_MODE_HALF         (CRASS_DEF_NUM_KMERS_4_MODE - (CRASS_DEF_NUM_KMERS_4_MODE/2)) // Ceil of 50% of CRASS_DEF_NUM_KMERS_4_MODE
#define CRASS_DEF_MIN_READ_DEPTH                2                   // read depth used for consensus building
#define CRASS_DEF_ZONE_EXT_CONS_CUT_OFF         (0.55)              // minimum identity to extend a DR from the "zone" outwards
#define CRASS_DEF_COLLAPSED_CONS_CUT_OFF        (0.75)              // minimum identity to identify a potential collapsed cluster
#define CRASS_DEF_COLLAPSED_THRESHOLD           (0.30)              // in the event that clustering has collapsed two DRs into one, this number is used to plait them apart
#define CRASS_DEF_PARTIAL_SIM_CUT_OFF           (0.85)              // The similarity needed to exted into partial matches
#define CRASS_DEF_DR_ZONE_TRIM_AMOUNT           (3)                 // trim this many bases off the start and end of the DR
#define CRASS_DEF_MIN_PARTIAL_LENGTH            (4)                 // The mininum length allowed for a partial direct repeat at the beginning or end of a read 
// --------------------------------------------------------------------
// GENOME ALGORITHM DEFINES
// --------------------------------------------------------------------
#define CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY  (0.82)
#define CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF     12
#define CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF     30
#define CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH         6
#define CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH         9
#define CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH     8
#define CRASS_DEF_SCAN_LENGTH                      30
#define CRASS_DEF_SCAN_CONFIDENCE                  (0.70)
#define CRASS_DEF_TRIM_EXTEND_CONFIDENCE           (0.5)
#define CRASS_DEF_DEFAULT_MIN_NUM_REPEATS          3
// --------------------------------------------------------------------
  // FILE IO
// --------------------------------------------------------------------
#define CRASS_DEF_FASTQ_FILENAME_MAX_LENGTH     1024
#define CRASS_DEF_DEF_KMER_LOOKUP_EXT           "crass_kmers.txt"
#define CRASS_DEF_DEF_PATTERN_LOOKUP_EXT        "crass_direct_repeats.txt"
#define CRASS_DEF_DEF_SPACER_LOOKUP_EXT         "crass_spacers.txt"
#define CRASS_DEF_CRISPR_EXT					".crispr"
// --------------------------------------------------------------------
// XML
// --------------------------------------------------------------------
#define CRASS_DEF_CRISPR_HEADER					"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<crass_assem version=\"1.0\">\n"
#define CRASS_DEF_ROOT_ELEMENT                  "crass_assem"
#define CRASS_DEF_XML_VERSION                   "1.0"
#define CRASS_DEF_CRISPR_FOOTER					"</crass_assem>\n"
// --------------------------------------------------------------------
// GRAPH BUILDING
// --------------------------------------------------------------------
#define CRASS_DEF_NODE_KMER_SIZE                7                  // size of the kmer that defines a crispr node
#define CRASS_DEF_MAX_CLEANING                  2                   // the maximum length that a branch can be befor it's cleaned
// --------------------------------------------------------------------
 // USER OPTION STRUCTURE
// --------------------------------------------------------------------
#define CRASS_DEF_STATS_REPORT                  false               // should we create a stats report
#define CRASS_DEF_STATS_REPORT_DELIM            "\t"                // delimiter string for stats report
#define CRASS_DEF_OUTPUT_DIR                    "./"                // default output directory
#define CRASS_DEF_MIN_DR_SIZE                   23                  // minimum direct repeat size
#define CRASS_DEF_MAX_DR_SIZE                   45                  // maximum direct repeat size
#define CRASS_DEF_MIN_SPACER_SIZE               26                  // minimum spacer size
#define CRASS_DEF_MAX_SPACER_SIZE               50                  // maximum spacer size
#define CRASS_DEF_NUM_DR_ERRORS                 0                   // maxiumum allowable errors in direct repeat
#define CRASS_DEF_COVCUTOFF                     3                   // minimum number of attached spacers that a group needs to have
#ifdef DEBUG
    #define CRASS_DEF_MAX_LOGGING               10
#else
    #define CRASS_DEF_MAX_LOGGING               4
#endif

#define CRASS_DEF_DEFAULT_LOGGING               1
#define CRASS_DEF_LOGTOSCREEN                   false               // should we log to screen or to file
#define CRASS_DEF_REMOVE_HOMOPOLYMERS           false               // should we remove homopolymers
#define CRASS_DEF_NO_SCALLING                   false               // perform scalling by default
#define CRASS_DEF_HOMOPOLYMER_SCALLING          (0.70)              // the scalling for the spacers and direct repeats when remove homopolymers is set
#define CRASS_DEF_NUM_OF_BINS                   -1                  // the number of bins to create
#define CRASS_DEF_GRAPH_COLOUR                  RED_BLUE            // default colour scale for the graphs
#define CRASS_DEF_SPACER_LONG_DESC              false               // use a long desc of the spacer in the output graph

typedef struct {
    int                 logLevel;                                           // level of verbosity allowed in the log file
    bool                reportStats;                                        // print a starts report currently not used
    unsigned int        lowDRsize;                                          // the lower size limit for a direct repeat
    unsigned int        highDRsize;                                         // the upper size limit for a direct repeat
    unsigned int        lowSpacerSize;                                      // the lower limit for a spacer
    unsigned int        highSpacerSize;                                     // the upper size limit for a spacer
    std::string         output_fastq;                                       // the output directory for the output files
    std::string         delim;                                              // delimiter used in stats report currently not used
    int                 kmer_size;                                          // number of kmers needed to be shared to add to a cluser
    unsigned int        searchWindowLength;                                 // option 'w'used in long read search only
    unsigned int        minNumRepeats;                                      // option 'n'used in long read search only
    bool                logToScreen;                                        // log to std::cout rather than to the log file
    bool                removeHomopolymers;                                 // correct for homopolymer errors
    int                 coverageBins;                                       // The number of bins of colours
    RB_TYPE             graphColourType;                                    // the colour type of the graph
    double              averageSpacerScalling;                              // decimal for reduction in the spacer size
    double              averageDrScalling;                                  // decimal for the reduction in the direct repeat size
    bool                dontPerformScalling;                                // turn all scalling off for the user to define variables
    std::string         layoutAlgorithm;                                    // the graphviz layout algorithm to use 
    bool                longDescription;                                    // print a long description for the final spacer graph
    int                 covCutoff;                                          // The lower bounds of acceptable numbers of reads that a group can have
} options;

// --------------------------------------------------------------------

#endif // __CRASSDEFINES_H
