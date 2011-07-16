/*
 *  crass.cpp is part of the CRisprASSembler project
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

// system includes
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string>
#include <unistd.h>
#include <fcntl.h>
#include <istream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>

// local includes
#include "crass.h"
#include "crass_defines.h"
#include "LoggerSimp.h"
#include "SeqUtils.h"
#include "WorkHorse.h"


//**************************************
// user input + system
//**************************************

void help_message() {
    std::cout<<"Usage:  "<<CRASS_DEF_PRG_NAME<<" [options]  <fastq_or_fasta_files>"<<std::endl<<std::endl;
    std::cout<< "\t-d --minDR <INT>           Minimim length of the direct repeat"<<std::endl; 
    std::cout<< "                             to search for [Default: 23]"<<std::endl;
    std::cout<< "\t-D --maxDR <INT>           Maximim length of the direct repeat"<<std::endl; 
    std::cout<<"                              to search for [Default: 47]"<<std::endl;
    std::cout<< "\t-h                         This help message"<<std::endl;
    std::cout<< "\t-k --kmerCount <INT>       The number of the kmers that need to be"<<std::endl; 
    std::cout<< "\t                           shared for clustering [Default: 8]"<<std::endl;
    std::cout<< "\t-l --logLevel <INT>        Output a log file and set a log level [1 - 10]"<<std::endl;
    std::cout<< "\t-m --maxMismatches <INT>   Total number of mismatches to at most allow for"<<std::endl;
    std::cout<< "\t                           in search pattern [Default: 0]"<<std::endl;
    std::cout<< "\t-o --outDir <DIRECTORY>    Output directory [default: .]"<<std::endl;
    std::cout<< "\t-s --minSpacer <INT>       Minimim length of the spacer to search for [Default: 26]"<<std::endl;
    std::cout<< "\t-S --maxSpacer <INT>       Maximim length of the spacer to search for [Default: 50]"<<std::endl;
    std::cout<< "\t-V --version               Program and version information"<<std::endl;
    std::cout<< "\t--dumpReads                Print reads containing CRISPR to file after the find stage"<<std::endl;

}

void version_info() {
    std::cout<<std::endl<<CRASS_DEF_LONG_NAME<<" ("<<CRASS_DEF_PRG_NAME<<")"<<std::endl<<"revison "<<CRASS_DEF_REVISION<<" version "<<CRASS_DEF_MAJOR_VERSION<<" subversion "<<CRASS_DEF_MINOR_VERSION<<" ("<<CRASS_DEF_VERSION<<")"<<std::endl<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Copyright (C) 2011 Connor Skennerton & Michael Imelfort"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY"<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it"<<std::endl;
    std::cout<<"under certain conditions: See the source for more details"<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
}

int process_options(int argc, char *argv[], options *opts) {
    int c;
    int index;
    //char *opt_o_value = NULL;
    char *opt_b_value = NULL;
    
    while( (c = getopt_long(argc, argv, "hVrl:k:p:m:o:b:cd:D:s:S:", long_options, &index)) != -1 ) {
        switch(c) {
            case 'h': help_message(); exit(0); break;
            case 'V': version_info(); exit(0); break;
            case 'o': opts->output_fastq = optarg; break;
            case 'p': opts->pat_file = optarg; break;
            case 'b': opt_b_value = optarg; break;
            case 'r': opts->report_stats = true; break;
            case 'm': opts->max_mismatches = atoi(optarg); break;
            case 'd': opts->lowDRsize = atoi(optarg); break;
            case 'D': opts->highDRsize = atoi(optarg); break;
            case 's': opts->lowSpacerSize = atoi(optarg); break;
            case 'S': opts->highSpacerSize = atoi(optarg); break;
            case 'c': opts->count = 1; break;
            case 'k': opts->kmer_size = atoi(optarg); break;
            case 'l': opts->logger_level = atoi(optarg); break;
            case '?': version_info(); help_message(); exit(1); break;
            case '0':
                if( strcmp( "dumpReads", long_options[index].name ) == 0 ) opts->detect = true;
                else if ( strcmp( "454", long_options[index].name ) == 0 ) opts->fourFiveFour = true;
                else if ( strcmp( "illumina", long_options[index].name ) == 0 ) opts->illumina = true;
                else if ( strcmp( "sanger", long_options[index].name ) == 0 ) opts->sanger = true;
                else if ( strcmp( "genome", long_options[index].name ) == 0 ) opts->genome = true;
                break;
            default: abort(); break;
        }
    }
    if (optind == argc)
    {
        std::cerr<<CRASS_DEF_PRG_NAME<<" : [ERROR] No input files were provided. Try ./"<<CRASS_DEF_PRG_NAME" -h for help."<<std::endl;
        exit(1);
        //        printf("no files given\n");
        //        exit(1);
    }
    
    /* setup delimiter for stats report (if given) */
    if ( opt_b_value != NULL ) 
    {
        strncpy(opts->delim, opt_b_value, CRASS_DEF_MAX_DELIM_LENGTH);
    }
    return optind;
}

//**************************************
// rock and roll
//**************************************

int main(int argc, char *argv[]) 
{
    
    int opt_idx;
    
    /* application of default options */
    options opts = {
        0,             // count flag
        false,         // exit after first find
        1,             // logging minimum by default
        false,             // illumina reads
        false,             // 454 reads
        false,             // sanger reads
        false,             // genome
        false,             // output stats report
        23,            // minimum direct repeat size
        47,            // maximum direct repeat size
        26,            // minimum spacer size
        50,            // maximum spacer size
        0,             // maxiumum allowable errors in direct repeat
        "./",            // output file directory
        "\t",          // delimiter string for stats report
        NULL,          //  pattern file name
        8             // the number of the kmers that need to be shared for clustering
    };
    
    // initialize the log file
    opt_idx = process_options(argc, argv, &opts);
    
    
    /* setup the appropriate output file pointer */
//    if ( opts.output_fastq.length() == 0 ) 
//    {
//        out_fp = "./";
//    }
//    else 
//    {
//        out_fp = opts.output_fastq;
//    }
    std::string logFileName = opts.output_fastq+ "crass.log";

    intialiseGlobalLogger(logFileName, opts.logger_level);


    if (opt_idx >= argc) 
    {
        std::cerr<<CRASS_DEF_PRG_NAME<<" : [ERROR] specify FASTQ files to process!"<<std::endl;
        exit(1);
    }
        
    // The remaining command line arguments are FASTQ(s) to process
    // put them in a vector to pass to the workhorse
    std::vector<std::string> seq_files;
    
    while (opt_idx < argc) 
    {
        //
        std::string seq_file(argv[opt_idx]);
        seq_files.push_back(seq_file);
        logTimeStamp();
        opt_idx++;
    }
    
    // Time to make the workhorse
    WorkHorse * mHorse = new WorkHorse(&opts);
    mHorse->doWork(seq_files);
    delete mHorse;
        
    return 0;
}
