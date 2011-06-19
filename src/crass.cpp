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

using namespace std;

//**************************************
// user input + system
//**************************************

void help_message() {
    fprintf(stdout, "Usage: %s %s %s\n", 
                    PRG_NAME, "[options]", "<fastq_or_fasta_files>");
    fprintf(stdout, "\t%-20s%-20s\n", "-d <INT>", "Minimim length of the direct repeat to search for [Default: 23]");
    fprintf(stdout, "\t%-20s%-20s\n", "-D <INT>", "Maximim length of the direct repeat to search for [Default: 47]");
    fprintf(stdout, "\t%-20s%-20s\n", "-h", "This help message");
    fprintf(stdout, "\t%-20s%-20s\n", "-k <INT>", "Length of the kmer used for assembling and indexing the");
    fprintf(stdout, "\t%-20s%-20s\n", "", "spacers and the direct repeats [Default: 12]");
    fprintf(stdout, "\t%-20s%-20s\n", "-l <INT>", "output a log file and set a log level [1 - 10]");
    fprintf(stdout, "\t%-20s%-20s\n", "-m <INT>", "Total number of mismatches to at most allow for");
    fprintf(stdout, "\t%-20s%-20s\n", "", "in search pattern [Default: 0]");
    fprintf(stdout, "\t%-20s%-20s\n", "-o <DIRECTORY>", "output directory [default: .]");
    fprintf(stdout, "\t%-20s%-20s\n", "-s <INT>", "Minimim length of the spacer to search for [Default: 26]");
    fprintf(stdout, "\t%-20s%-20s\n", "-S <INT>", "Maximim length of the spacer to search for [Default: 50]");
    fprintf(stdout, "\t%-20s%-20s\n", "-V", "Program and version information");

}

void version_info() {
    cout<<endl<<LONG_NAME<<" ("<<PRG_NAME<<")"<<endl<<"revison "<<MCD_REVISION<<" version "<<MCD_MAJOR_VERSION<<" subversion "<<MCD_MINOR_VERSION<<" ("<<MCD_VERSION<<")"<<endl<<endl;
    cout<<"---------------------------------------------------------------"<<endl;
    cout<<"Copyright (C) 2011 Connor Skennerton & Michael Imelfort"<<endl;
    cout<<"This program comes with ABSOLUTELY NO WARRANTY"<<endl;
    cout<<"This is free software, and you are welcome to redistribute it"<<endl;
    cout<<"under certain conditions: See the source for more details"<<endl;
    cout<<"---------------------------------------------------------------"<<endl;
}

int process_options(int argc, char *argv[], options *opts) {
    int c;
    //int index;
    //char *opt_o_value = NULL;
    char *opt_b_value = NULL;
    
    while( (c = getopt(argc, argv, "hVfrl:evak:p:m:o:b:cd:D:s:S:")) != -1 ) {
        switch(c) {
            case 'h':
                help_message();
                exit(0);
                break;
            case 'V':
                version_info();
                exit(0);
                break;
            case 'v':
                opts->invert_match = 1;
                break;
            case 'a':
                opts->show_all_records = 1;
                break;
            case 'o':
                opts->output_fastq = optarg;
                break;
            case 'p':
                opts->pat_file = optarg;
                break;
            case 'b':
                opt_b_value = optarg;
                break;
            case 'f':
                opts->report_fasta = 1;
                opts->report_fastq = 0;
                opts->report_stats = 0;
                break;
            case 'r':
                opts->report_fasta = 0;
                opts->report_fastq = 0;
                opts->report_stats = 1;
                break;
            case 'm':
                opts->max_mismatches = atoi(optarg);
                break;
            case 'd':
                opts->lowDRsize = atoi(optarg);
                break;
            case 'D':
                opts->highDRsize = atoi(optarg);
                break;
            case 's':
                opts->lowSpacerSize = atoi(optarg);
                break;
            case 'S':
                opts->highSpacerSize = atoi(optarg);
                break;
            case 'c':
                opts->count = 1;
                break;
            case 'e':
                opts->detect = true;
                break;
            case 'k':
                opts->kmer_size = atoi(optarg);
                break;
            case 'l':
                opts->logger_level = atoi(optarg);
                break;
            case '?':
                exit(1);
            default:
                abort();
        }
    }
    if (optind == argc)
    {
        fprintf(stderr, "%s : [ERROR] No input files were provided. Try ./%s -h for help.\n", PRG_NAME, PRG_NAME);
        exit(1);
        //        printf("no files given\n");
        //        exit(1);
    }
    
    /* setup delimiter for stats report (if given) */
    if ( opt_b_value != NULL ) 
    {
        strncpy(opts->delim, opt_b_value, DEF_MCD_MAX_DELIM_LENGTH);
    }
    return optind;
}

//**************************************
// rock and roll
//**************************************

int main(int argc, char *argv[]) 
{
    
    int opt_idx;
    std::string out_fp;
    
    /* application of default options */
    options opts = {
        0,             // count flag
        false,         // exit after first find
        1,             // logging off by default
        0,             // invert match flag
        0,             // show all records flag
        1,             // output fastq report
        0,             // output fasta report
        0,             // output stats report
        23,            // minimum direct repeat size
        47,            // maximum direct repeat size
        26,            // minimum spacer size
        50,            // maximum spacer size
        0,             // maxiumum allowable errors in direct repeat
        "",            // output file name
        "\t",          // delimiter string for stats report
        NULL,          //  pattern file name
        12             // the length of the kmer for ordering and assembling
    };
    
    // initialize the log file
    opt_idx = process_options(argc, argv, &opts);
    
    
    /* setup the appropriate output file pointer */
    if ( opts.output_fastq.length() == 0 ) 
    {
        out_fp = "./";
    }
    else 
    {
        out_fp = opts.output_fastq;
    }
    string logFileName = out_fp+ "mcd.log";

    intialiseGlobalLogger(logFileName, opts.logger_level);


    if (opt_idx >= argc) 
    {
        fprintf(stderr, "%s : %s\n",
                PRG_NAME, "[ERROR] specify FASTQ files to process!");
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
    WorkHorse * mHorse = new WorkHorse(&opts, out_fp);
    mHorse->doWork(seq_files);
    delete mHorse;
        
    return 0;
}
