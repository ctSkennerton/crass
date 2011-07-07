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
    std::string output_fastq;
    char delim[CRASS_DEF_MAX_DELIM_LENGTH];                         // delimiter used in stats report
    char * pat_file;
    int kmer_size;
} options;
// --------------------------------------------------------------------
#endif // __CRASSDEFINES_H
