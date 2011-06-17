//
//  MCD_defines.h
//  MCD
//
//  Created by Connor Skennerton on 7/05/11.
//  Copyright 2011 The Faculty of EAIT. All rights reserved.
//
#ifndef __MCDDEFINES_H
    #define __MCDDEFINES_H

#include <string>

// used in all parts of the program
#define LONG_NAME "MetagenomicCRISPRDetector"
#define MCD_VERSION "0.0.1"
#define MCD_MAJOR_VERSION 0
#define MCD_MINOR_VERSION 1
#define MCD_REVISION 0
#define PRG_NAME "mcd"

#define DEF_MCD_LOOKUP_EXT ".txt"
#define DEF_MCD_SPACER_LOOKUP "_spacers"
#define DEF_MCD_KMERS_LOOKUP "_kmers"
#define DEF_MCD_PATTERNS_LOOKUP "_direct_repeats"
#define DEF_MCD_SPACER_KMERS_LOOKUP "_spacer_kmers"

// the length of read where the multipattern matcher would become insignificant
// ie where it is very likely that there is 2 full repeat-spacer-repeat units
#define DEF_MCD_READ_LENGTH_CUTOFF 200.0

#define DEF_MCD_FASTQ_FILENAME_MAX_LENGTH 1024
#define DEF_MCD_MAX_PATTERN_LENGTH 1024
#define DEF_MCD_MAX_DELIM_LENGTH 10

// user options structure
typedef struct {
    int count;
    bool detect;
    int logger_level;
    int invert_match;
    int show_all_records;
    int report_fastq;
    int report_fasta;
    int report_stats;
    int lowDRsize;
    int highDRsize;
    int lowSpacerSize;
    int highSpacerSize;
    int max_mismatches;
    //    int max_deletions;
    //    int max_substitutions;
    std::string output_fastq;
    //    char search_pattern[DEF_MCD_MAX_PATTERN_LENGTH];
    char delim[DEF_MCD_MAX_DELIM_LENGTH];         /* delimiter used in stats report */
    char * pat_file;
    int kmer_size;
} options;

// File IO
#define MCD_DEF_KMER_LOOKUP_EXT         "mcd_kmers.txt"
#define MCD_DEF_PATTERN_LOOKUP_EXT      "mcd_direct_repeats.txt"
#define MCD_DEF_SPACER_LOOKUP_EXT       "mcd_spacers.txt"


#endif
