
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
#include <time.h>

// local includes
#include <config.h>
#include "crass.h"
#include "crass_defines.h"
#include "LoggerSimp.h"
#include "SeqUtils.h"
#include "WorkHorse.h"

// added for testing only
#include "NodeManager.h"
#include "SpacerInstance.h"
#include "CrisprNode.h"
#include "ReadHolder.h"
#include "StringCheck.h"


//**************************************
// user input + system
//**************************************

void help(void)
{
    std::cout<<"Usage:  "<<PACKAGE_NAME<<"  [options] { <inputFiles> <read_type> } {...}"<<std::endl<<std::endl;
    std::cout<<"for further information try:"<<std::endl;
    std::cout<<PACKAGE_NAME<<" -h"<<std::endl;
}



void usage(void) 
{
    std::cout<<"Usage:  "<<PACKAGE_NAME<<"  [options] { <inputFiles> <read_type> } {...}"<<std::endl<<std::endl;
    std::cout<<"where read_type is one of the following:"<<std::endl;
    std::cout<<"--sanger"<<std::endl<<"--illumina"<<std::endl<<"--454"<<std::endl<<"--ionTorrent"<<std::endl<<std::endl;
    std::cout<<"Options:"<<std::endl;
    std::cout<< "\t-d --minDR <INT>           Minimim length of the direct repeat"<<std::endl; 
    std::cout<< "\t                           to search for [Default: 23]"<<std::endl;
    std::cout<< "\t-D --maxDR <INT>           Maximim length of the direct repeat"<<std::endl; 
    std::cout<< "\t                            to search for [Default: 47]"<<std::endl;
    std::cout<< "\t-h --help                  This help message"<<std::endl;
    std::cout<< "\t-k --kmerCount <INT>       The number of the kmers that need to be"<<std::endl; 
    std::cout<< "\t                           shared for clustering [Default: 8]"<<std::endl;
//#ifdef DEBUG
    std::cout<< "\t-l --logLevel <INT>        Output a log file and set a log level [1 - 10]"<<std::endl;
//#else
//    std::cout<< "\t-l --logLevel <INT>        Output a log file and set a log level [1 - 4]"<<std::endl;
//#endif
//    std::cout<< "\t-m --maxMismatches <INT>   Total number of mismatches to at most allow for"<<std::endl;
//    std::cout<< "\t                           in search pattern [Default: 0]"<<std::endl;
    std::cout<< "\t-n --minNumRepeats <INT>   Total number of direct repeats in a CRISPR for"<<std::endl;
    std::cout<< "\t                           it to be considered real [Default: 3]"<<std::endl;
    std::cout<< "\t-o --outDir <DIRECTORY>    Output directory [default: .]"<<std::endl;
    std::cout<< "\t-s --minSpacer <INT>       Minimim length of the spacer to search for [Default: 26]"<<std::endl;
    std::cout<< "\t-S --maxSpacer <INT>       Maximim length of the spacer to search for [Default: 50]"<<std::endl;
    std::cout<< "\t-V --version               Program and version information"<<std::endl;
    std::cout<< "\t-w --windowLength <INT>    The length of the search window. Can only be"<<std::endl; 
    std::cout<< "\t                           a number between 6 - 9 [Default: 8]"<<std::endl;
    std::cout<< "\t--dumpReads                Print reads containing CRISPR to file after the find stage"<<std::endl;
    std::cout<< "\t--logToScreen              Print the logging information to screen rather than a file"<<std::endl;
}

void versionInfo(void) 
{
    std::cout<<std::endl<<PACKAGE_FULL_NAME<<" ("<<PACKAGE_NAME<<")"<<std::endl<<"version "<<PACKAGE_MAJOR_VERSION<<" subversion "<<PACKAGE_MINOR_VERSION<<" revison "<<PACKAGE_REVISION<<" ("<<PACKAGE_VERSION<<")"<<std::endl<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Copyright (C) 2011 Connor Skennerton & Michael Imelfort"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY"<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it"<<std::endl;
    std::cout<<"under certain conditions: See the source for more details"<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
}

int processOptions(int argc, char *argv[], options *opts) 
{
    int c;
    int index;
    int l_val;
    int w_val;
    //char *opt_o_value = NULL;
    char *opt_b_value = NULL;
    
    while( (c = getopt_long(argc, argv, "+w:hVrl:k:p:n:o:b:cd:D:s:S:", long_options, &index)) != -1 ) {
        switch(c) {
            case 'n': opts->minNumRepeats = atoi(optarg); break;
            case 'w': 
                w_val = atoi(optarg); 
                if ((w_val < CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH) || (w_val > CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH))
                {
                    // Change window length
                    opts->searchWindowLength = CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH;
                    std::cerr<<"WARNING: Specified window length higher than max. Changing window length to " << opts->searchWindowLength << " instead of " << w_val<<std::endl;
                }
                break;        
            case 'h': versionInfo(); usage();exit(0); break;
            case 'V': versionInfo(); exit(0); break;
            case 'o': opts->output_fastq = optarg; break;
            case 'p': opts->pat_file = optarg; break;
            case 'b': opt_b_value = optarg; break;
            case 'r': opts->report_stats = true; break;
//            case 'm': opts->max_mismatches = atoi(optarg); break;
            case 'd': opts->lowDRsize = atoi(optarg); break;
            case 'D': opts->highDRsize = atoi(optarg); break;
            case 's': opts->lowSpacerSize = atoi(optarg); break;
            case 'S': opts->highSpacerSize = atoi(optarg); break;
            case 'c': opts->count = 1; break;
            case 'k': opts->kmer_size = atoi(optarg); break;
            case 'l': 
                l_val =  atoi(optarg);
//#ifdef DEBUG
                if (l_val > CRASS_DEF_MAX_DEBUG_LOGGING)
                {
                    std::cerr<<"WARNING: Specified log level higher than max. Changing log level to "<<CRASS_DEF_MAX_DEBUG_LOGGING<<" instead of "<<l_val<<std::endl;
                    opts->logger_level = CRASS_DEF_MAX_DEBUG_LOGGING;
                }
                else
                {
                    opts->logger_level = l_val;
                }
//#else
//                if(l_val > CRASS_DEF_MAX_LOGGING)
//                {
//                    std::cerr<<"WARNING: Specified log level higher than max. Changing log level to "<<CRASS_DEF_MAX_LOGGING<<" instead of "<<l_val<<std::endl;
//                    opts->logger_level = CRASS_DEF_MAX_LOGGING;
//                }
//                else
//                {
//                    opts->logger_level = l_val;
//                }
//#endif
                break;
            case 0:
                if( strcmp( "dumpReads", long_options[index].name ) == 0 ) opts->detect = true;
//                else if ( strcmp( "454", long_options[index].name ) == 0 ) opts->fourFiveFour = true;
//                else if ( strcmp( "illumina", long_options[index].name ) == 0 ) opts->illumina = true;
//                else if ( strcmp( "sanger", long_options[index].name ) == 0 ) opts->sanger = true;
//                else if ( strcmp( "ionTorrent", long_options[index].name ) == 0 ) opts->fourFiveFour = true;
                else if ( strcmp( "logToScreen", long_options[index].name ) == 0 ) opts->logToScreen = true;

                break;
            default: versionInfo(); help(); exit(1); break;
        }
    }
//    if (optind == argc)
//    {
//        std::cerr<<PACKAGE_NAME<<" : [ERROR] No input files were provided. Try ./"<<PACKAGE_NAME" -h for help."<<std::endl;
//        exit(1);
//    }
    
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
    
    if(argc == 1) 
    {
        help();
        exit(1);
    }
    /* application of default options */
    options opts = {
        0,                                      // count flag
        false,                                  // exit after first find
        CRASS_DEF_DEFAULT_LOGGING,              // logging minimum by default
        false,                                  // output stats report
        CRASS_DEF_MIN_DR_SIZE,                  // minimum direct repeat size
        CRASS_DEF_MAX_DR_SIZE,                  // maximum direct repeat size
        CRASS_DEF_MIN_SPACER_SIZE,              // minimum spacer size
        CRASS_DEF_MAX_SPACER_SIZE,              // maximum spacer size
        "./",                                   // output file directory
        "\t",                                   // delimiter string for stats report
        NULL,                                   //  pattern file name
        CRASS_DEF_K_CLUST_MIN,                  // the number of the kmers that need to be shared for clustering
        CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH, // the search window length
        CRASS_DEF_DEFAULT_MIN_NUM_REPEATS,      // mininum number of repeats for long read algorithm 
        false                                   // should we log to screen
    };

    int opt_idx = processOptions(argc, argv, &opts);

    if (opt_idx >= argc) 
    {
        std::cerr<<PACKAGE_NAME<<" : [ERROR] Specify sequence files to process!"<<std::endl;
        help();
        exit(1);
    }
    
    // initialize the log file
    std::string logFileName;
    if (opts.logToScreen) 
    {
        logFileName = "";
    } 
    else 
    {
        // create a time stamp for the log file
        time_t rawtime;
        struct tm * timeinfo;
        char buffer [80];
        
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        
        strftime (buffer,80,"%d_%m_%Y_%H%M%S",timeinfo);
        std::string tmp(buffer);
        logFileName = opts.output_fastq+ "crass."+tmp+".log";
        
    }
    
    intialiseGlobalLogger(logFileName, opts.logger_level);
    
    // The remaining command line arguments are FASTQ(s) to process
    // put them in a vector to pass to the workhorse
    std::vector<std::string> seq_files;
    std::vector<SequenceType> seq_type_files;
    
    while (opt_idx < argc) 
    {
        //
        std::string seq_file(argv[opt_idx]);
        
        seq_files.push_back(seq_file);
        
        opt_idx++;

        if (opt_idx >= argc || '-' != argv[opt_idx][0]) 
        {
            //throw error
            std::cerr<<PACKAGE_NAME<<" : [ERROR] A file must be followed by a read type!"<<std::endl;
            help();
            exit(1);
            
        } 
        else 
        {
            //check for the sequence type
            if (strcmp(argv[opt_idx],"--454") == 0) 
            {
                seq_type_files.push_back(fourFiveFour);
            }
            else if(strcmp(argv[opt_idx],"--illumina") == 0)
            {
                seq_type_files.push_back(illumina);
            }
            else if(strcmp(argv[opt_idx],"--sanger") == 0)
            {
                seq_type_files.push_back(sanger);
            }
            else if(strcmp(argv[opt_idx],"--ionTorrent") == 0)
            {
                seq_type_files.push_back(fourFiveFour);
            }
            else
            {
                // throw error
                std::cerr<<PACKAGE_NAME<<" : [ERROR] The sequence format you provided, "<<argv[opt_idx]<<", is not recognized!"<<std::endl;
                help();
                exit(1);
            }
        }
        opt_idx++;
    }
    logTimeStamp();
    
    std::string test_read = "CCGTGCTACCCGCAGGACAATACGGCCCCATGGCCGCGTTCCCCGCAGGCGCGGGGATGAACCGCGGGCGCGGGTGCGGGAGGCTTTGCGACGCAGGCGTTCCCCGCAGGCGCGGGGATGAACCGAGCACGGTCAACGCAGGAGGCACCGAGTACCGGCGTTCCCCGCAGGCGCGGGGATGAACCGTCTTTATCGGTAGCTGAGGCTCTTGCGCAAAGGCGTTCCCCGCAGGCGCGGGGATGAACCGAGCATTCAGGCTAAACACTGAGCACTAATGATGCGTTCCCCGCAGGCGCGGGGATGAACCGCCGGAGTAGAGGCAGAAGACGAGACGCCACAAGCGTTCCCCGCAGGCGCGGGGATGAACCGGGTGATAGAGCTT";
    std::string test_header = "test_1";
    ReadHolder * rh = new ReadHolder(test_read,test_header);
    rh->add(35, 63);
    rh->add(96, 124);
    rh->add(157, 185);
    rh->add(218, 246);
    rh->add(279, 307);
    
    StringCheck * sc = new StringCheck();
    
    NodeManager * nm = new NodeManager("GCGTTCCCCGCAGGCGCGGGGATGAACCG", sc);
    
    nm->addReadHolder(rh);
    
    NodeListIterator nm_iter = nm->nodeBegin();
    while (nm_iter != nm->nodeEnd()) 
    {
        (*nm_iter)->printEdges();
        nm_iter++;
    }
    
    return 0;
    
//    WorkHorse * mHorse = new WorkHorse(&opts);
//    mHorse->doWork(seq_files, seq_type_files);
//    delete mHorse;
//    return 0;
}
