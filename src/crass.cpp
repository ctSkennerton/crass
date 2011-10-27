
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
#include <sys/types.h>
#include <sys/stat.h>

// local includes
#include "config.h"
#include "crass.h"
#include "crassDefines.h"
#include "LoggerSimp.h"
#include "SeqUtils.h"
#include "WorkHorse.h"
#include "Rainbow.h"
#include "AssemblyWrapper.h"


//**************************************
// user input + system
//**************************************

void usage(void) 
{
    std::cout<<"Usage:  "<<PACKAGE_NAME<<"  [options] { inputFile ...}"<<std::endl<<std::endl;
    std::cout<<"General Options:"<<std::endl;
    std::cout<< "-h --help                    This help message"<<std::endl;
    std::cout<< "-l --logLevel        <INT>   Output a log file and set a log level [1 - "<<CRASS_DEF_MAX_LOGGING<<"]"<<std::endl;
    std::cout<< "-o --outDir          <DIR>   Output directory [default: .]"<<std::endl;
    std::cout<< "-V --version                 Program and version information"<<std::endl;
    std::cout<< "--logToScreen                Print the logging information to screen rather than a file"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"CRISPR Identification Options:"<<std::endl;
    std::cout<< "-d --minDR           <INT>   Minimim length of the direct repeat"<<std::endl; 
    std::cout<< "                             to search for [Default: "<<CRASS_DEF_MIN_DR_SIZE<<"]"<<std::endl;
    std::cout<< "-D --maxDR           <INT>   Maximim length of the direct repeat"<<std::endl; 
    std::cout<< "                             to search for [Default: "<<CRASS_DEF_MAX_DR_SIZE<<"]"<<std::endl;
    std::cout<< "-n --minNumRepeats   <INT>   Total number of direct repeats in a CRISPR for"<<std::endl;
    std::cout<< "                             it to be considered real [Default: "<<CRASS_DEF_DEFAULT_MIN_NUM_REPEATS<<"]"<<std::endl;
    std::cout<< "-s --minSpacer       <INT>   Minimim length of the spacer to search for [Default: "<<CRASS_DEF_MIN_SPACER_SIZE<<"]"<<std::endl;
    std::cout<< "-S --maxSpacer       <INT>   Maximim length of the spacer to search for [Default: "<<CRASS_DEF_MAX_SPACER_SIZE<<"]"<<std::endl;
    std::cout<< "-w --windowLength    <INT>   The length of the search window. Can only be"<<std::endl; 
    std::cout<< "                             a number between "<<CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH<<" - "<<CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH<<" [Default: "<<CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH<<"]"<<std::endl;
    std::cout<< "-x --spacerScalling  <REAL>  A decimal number that represents the reduction in size of the spacer"<<std::endl;
    std::cout<< "                             when the --removeHomopolymers option is set [Default: "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<"]"<<std::endl;
    std::cout<< "-y --repeatScalling  <REAL>  A decimal number that represents the reduction in size of the direct repeat"<<std::endl;
    std::cout<< "                             when the --removeHomopolymers option is set [Default: "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<"]"<<std::endl;
    std::cout<< "--noScalling                 Use the given spacer and direct repeat ranges when --removeHomopolymers is set. "<<std::endl;
    std::cout<< "                             The default is to scale the numbers by "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<" or by values set using -x or -y"<<std::endl;
    std::cout<< "--dumpReads                  Print reads containing CRISPR to file after the find stage"<<std::endl;
    std::cout<< "--removeHomopolymers         Correct for homopolymer errors [default: no correction]"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"CRISPR Assembly Options:"<<std::endl;
    std::cout<< "-k --kmerCount       <INT>   The number of the kmers that need to be"<<std::endl; 
    std::cout<< "                             shared for clustering [Default: "<<CRASS_DEF_K_CLUST_MIN<<"]"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Graph Output Options: "<<std::endl;
    std::cout<<"-b --numBins          <INT>   The number of colour bins for the output graph."<<std::endl;
    std::cout<<"                              Default is to have as many colours as there are"<<std::endl;
    std::cout<<"                              different values for the coverage of Nodes in the graph"<<std::endl;
    std::cout<<"-c --graphColour     <TYPE>   Defines the range of colours to use for the output graph"<<std::endl;
    std::cout<<"                              the different types available are:"<<std::endl;
    std::cout<<"                              red-blue, blue-red, green-red-blue, red-blue-green"<<std::endl;
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
    bool scalling = false;
    while( (c = getopt_long(argc, argv, "w:hVrl:k:n:o:b:c:d:D:s:S:x:y:", long_options, &index)) != -1 ) 
    {
        switch(c) 
        {
            case 'n':
                if (atoi(optarg) < 2) 
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: The mininum number of repeats cannot be less than 2"<<std::endl;
                    usage();
                    exit(1);
                }
                opts->minNumRepeats = atoi(optarg); 
                break;
            case 'w': 
                opts->searchWindowLength = atoi(optarg); 
                if ((opts->searchWindowLength < CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH) || (opts->searchWindowLength > CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH))
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified window length higher than max. Changing window length to " << CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH << " instead of " << opts->searchWindowLength<<std::endl;
                    // Change window length
                    opts->searchWindowLength = CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH;
                }
                break;        
            case 'h': 
                versionInfo(); 
                usage();
                exit(0); 
                break;
            case 'V': 
                versionInfo(); 
                exit(0); 
                break;
            case 'o': 
                opts->output_fastq = optarg; 
                // just in case the user put '.' or '..' or '~' as the output directory
                if (opts->output_fastq[opts->output_fastq.length() - 1] != '/')
                {
                    opts->output_fastq += '/';
                }
                
                // check if our output folder exists
                struct stat file_stats;
                if (stat(opts->output_fastq.c_str(),&file_stats)) 
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: The Directory "<<opts->output_fastq<<" does not exist"<<std::endl;
                    usage();
                    exit(1);
                }
                break;
            case 'b': 
                if (atoi(optarg) < 1) 
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: The number of bins for coverage cannot be less than 1"<<std::endl;
                    usage();
                    exit(1);
                }
                opts->coverageBins = atoi(optarg); 
                break;
            case 'c': 
                if (strcmp(optarg, "red-blue") == 0) 
                {
                    opts->graphColourType = RED_BLUE;
                } 
                else if(strcmp(optarg, "read-blue-green") == 0)
                {
                    opts->graphColourType = RED_BLUE_GREEN;
                }
                else if(strcmp(optarg, "blue-red") == 0)
                {
                    opts->graphColourType = BLUE_RED; 
                }
                else if(strcmp(optarg, "green-blue-red") == 0)
                {
                    opts->graphColourType = GREEN_BLUE_RED;
                }
                else
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Unknown graph colour type "<<optarg<<" changing to default colour type (red-blue)"<<std::endl;
                    opts->graphColourType = RED_BLUE;
                }
                break;
            case 'r': 
                opts->reportStats = true; 
                break;
            case 'd': 
                opts->lowDRsize = atoi(optarg); 
                if (opts->lowDRsize < 8) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for direct repeat sizes cannot be "<<opts->lowDRsize<<" changing to "<<CRASS_DEF_MIN_DR_SIZE<<std::endl;
                    opts->lowDRsize = CRASS_DEF_MIN_DR_SIZE;
                }
                break;
            case 'D': 
                opts->highDRsize = atoi(optarg); 
                break;
            case 's': 
                opts->lowSpacerSize = atoi(optarg); 
                if (opts->lowSpacerSize < 8) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for spacer sizes cannot be "<<opts->lowSpacerSize<<" changing to "<<CRASS_DEF_MIN_SPACER_SIZE<<std::endl;
                    opts->lowSpacerSize = CRASS_DEF_MIN_SPACER_SIZE;
                }
                break;
            case 'S': 
                opts->highSpacerSize = atoi(optarg); 
                break;
            case 'k': 
                opts->kmer_size = atoi(optarg);
                if (opts->kmer_size < 6) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for kmer clustering cannot be "<<opts->kmer_size<<" changing to "<<CRASS_DEF_K_CLUST_MIN<<std::endl;
                    opts->kmer_size = CRASS_DEF_K_CLUST_MIN;
                }
                break;
            case 'l': 
                opts->logLevel =  atoi(optarg);
                if(opts->logLevel > CRASS_DEF_MAX_LOGGING)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified log level higher than max. Changing log level to "<<CRASS_DEF_MAX_LOGGING<<" instead of "<<opts->logLevel<<std::endl;
                    opts->logLevel = CRASS_DEF_MAX_LOGGING;
                }
                break;
            case 'x':
                opts->averageSpacerScalling = atof(optarg);
                if (isNotDecimal(opts->averageSpacerScalling)) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The average spacer scalling must be a decimal number. Changing to "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<" instead of "<<opts->averageSpacerScalling<<std::endl;
                    opts->averageSpacerScalling = CRASS_DEF_HOMOPOLYMER_SCALLING;
                }
                scalling = true;
                break;
            case 'y':
                opts->averageDrScalling = atof(optarg);
                if (isNotDecimal(opts->averageDrScalling)) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The average spacer scalling must be a decimal number. Changing to "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<" instead of "<<opts->averageDrScalling<<std::endl;
                    opts->averageDrScalling = CRASS_DEF_HOMOPOLYMER_SCALLING;
                }
                scalling = true;
                break;
            case 0:
                if( strcmp( "dumpReads", long_options[index].name ) == 0 ) opts->detect = true;
                else if ( strcmp( "removeHomopolymers", long_options[index].name ) == 0 ) opts->removeHomopolymers = true;
                else if ( strcmp( "logToScreen", long_options[index].name ) == 0 ) opts->logToScreen = true;
                else if ( strcmp( "noScalling", long_options[index].name ) == 0 ) opts->dontPerformScalling = true;
                break;
            default: 
                versionInfo(); 
                usage(); 
                exit(1); 
                break;
        }
    }
    // Sanity checks for the high and low dr size
    if (opts->lowDRsize >= opts->highDRsize) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: The lower direct repeat bound is bigger than the higher bound ("<<opts->lowDRsize<<" >= "<<opts->highDRsize<<")"<<std::endl;
        usage();
        exit(1);
    }
    // Sanity checks for the high and low spacer size
    if (opts->lowSpacerSize >= opts->highSpacerSize) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: The lower spacer bound is bigger than the higher bound ("<<opts->lowSpacerSize<<" >= "<<opts->highSpacerSize<<")"<<std::endl;
        usage();
        exit(1);
    }
    
    // sanity check so that the user doesn't specify scalling and no scalling simultaneously
    if (scalling & opts->dontPerformScalling) {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: Cannot use scalling (-x -y) in conjunction with --noScalling"<<std::endl;
        usage();
        exit(1);
    }
    
    // warn them if they try to scale without specifying to remove homopolymers
    if (scalling && !opts->removeHomopolymers) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: scalling (-x -y) can only be used in conjunction with --removeHomopolymers"<<std::endl;
        usage();
        exit(1);
    }
    
    // scale the direct repeat and spacer lengths if we should
    if (opts->removeHomopolymers && !opts->dontPerformScalling) 
    {
        opts->lowDRsize *= opts->averageDrScalling;
        opts->highDRsize *= opts->averageDrScalling;
        opts->lowSpacerSize *= opts->averageSpacerScalling;
        opts->highSpacerSize *= opts->averageSpacerScalling;
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
        usage();
        exit(1);
    }
    
    if (strcmp(argv[1], "assemble") == 0) {
        // our user wants to do an assembly so let's load up the assembler main function
        assemblyMain(argc - 1 , argv + 1);
        
    }
    /* application of default options */
    options opts = {
        false,                                  // exit after first find
        CRASS_DEF_DEFAULT_LOGGING,              // logging minimum by default
        false,                                  // output stats report
        CRASS_DEF_MIN_DR_SIZE,                  // minimum direct repeat size
        CRASS_DEF_MAX_DR_SIZE,                  // maximum direct repeat size
        CRASS_DEF_MIN_SPACER_SIZE,              // minimum spacer size
        CRASS_DEF_MAX_SPACER_SIZE,              // maximum spacer size
        "./",                                   // output file directory
        CRASS_DEF_STATS_REPORT_DELIM,           // delimiter string for stats report
        CRASS_DEF_K_CLUST_MIN,                  // the number of the kmers that need to be shared for clustering
        CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH, // the search window length
        CRASS_DEF_DEFAULT_MIN_NUM_REPEATS,      // mininum number of repeats for long read algorithm 
        false,                                  // should we log to screen
        false,                                  // should we remove homopolymers in reads
        -1,                                     // default for the number of bins of colours to create
        RED_BLUE,                               // default colour type of the graph
        CRASS_DEF_HOMOPOLYMER_SCALLING,         // average spacer scalling
        CRASS_DEF_HOMOPOLYMER_SCALLING,         // average direct repeat scalling
        false                                   // perform scalling by default
    };

    int opt_idx = processOptions(argc, argv, &opts);

    if (opt_idx >= argc) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: Specify sequence files to process!"<<std::endl;
        usage();
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
    
    intialiseGlobalLogger(logFileName, opts.logLevel);
    
    // The remaining command line arguments are FASTQ(s) to process
    // put them in a vector to pass to the workhorse
    std::vector<std::string> seq_files;
    
    while (opt_idx < argc) 
    {
        //
        std::string seq_file(argv[opt_idx]);
        
        seq_files.push_back(seq_file);
        
        opt_idx++;
    }
    logTimeStamp();
    
    WorkHorse * mHorse = new WorkHorse(&opts);
    mHorse->doWork(seq_files);
    delete mHorse;
    return 0;
}
