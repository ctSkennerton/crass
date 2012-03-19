
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
#include "WorkHorse.h"
#include "Rainbow.h"
#include "StlExt.h"
#ifdef PERFORM_CRASS_ASSEMBLY
    #include "AssemblyWrapper.h"
#endif


//**************************************
// user input + system
//**************************************

void usage(void) 
{
    std::cout<<"Compiler Options:"<<std::endl;
    std::cout<<"RENDERING =";
#ifdef RENDERING
    std::cout<<" 1"<<std::endl;
#else
    std::cout<<" 0"<<std::endl;
#endif
    std::cout<<"DEBUG =";
#ifdef DEBUG
    std::cout<<" 1"<<std::endl;
#else
    std::cout<<" 0"<<std::endl;
#endif
    std::cout<<"MEMCHECK =";
#ifdef MEMCHECK
    std::cout<<" 1"<<std::endl;
#else
    std::cout<<" 0"<<std::endl;
#endif
    std::cout<<"ASSEMBER =";
#ifdef PERFORM_CRASS_ASSEMBLY
    std::cout<<" 1"<<std::endl;
#else
    std::cout<<" 0"<<std::endl;
#endif
    std::cout<<"VERBOSE_LOGGER =";
#ifdef SUPER_LOGGING
    std::cout<<" 1"<<std::endl;
#else
    std::cout<<" 0"<<std::endl;
#endif
    std::cout<<std::endl;
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
    std::cout<< "--removeHomopolymers         Correct for homopolymer errors [default: no correction]"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"CRISPR Assembly Options:"<<std::endl;
    std::cout<< "-f --covCutoff       <INT>   Remove groups with less than x spacers [Default: "<<CRASS_DEF_COVCUTOFF<<"]"<<std::endl;
    std::cout<< "-k --kmerCount       <INT>   The number of the kmers that need to be"<<std::endl; 
    std::cout<< "                             shared for clustering [Default: "<<CRASS_DEF_K_CLUST_MIN<<"]"<<std::endl;
    std::cout<< "-K --graphNodeLen    <INT>   Length of the kmers used to make crispr nodes [Default: "<<CRASS_DEF_NODE_KMER_SIZE<<"]"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Output Options: "<<std::endl;
#ifdef RENDERING
    std::cout<<"-a --layoutAlgorithm  <TYPE>  Graphviz layout algorithm to use for printing spacer graphs. The following are available:"<<std::endl;
    #ifdef HAVE_NEATO
    std::cout<<"                              \tneato"<<std::endl;
    #endif
    #ifdef HAVE_DOT
    std::cout<<"                              \tdot [Default]"<<std::endl;
    #endif
    #ifdef HAVE_FDP
    std::cout<<"                              \tfdp"<<std::endl;
    #endif
    #ifdef HAVE_SFDP
    std::cout<<"                              \tsfdp"<<std::endl;
    #endif
    #ifdef HAVE_CIRCO
    std::cout<<"                              \tcirco"<<std::endl;
    #endif
    #ifdef HAVE_TWOPI
    std::cout<<"                              \ttwopi"<<std::endl;
    #endif
#endif 
    std::cout<<"-b --numBins          <INT>   The number of colour bins for the output graph."<<std::endl;
    std::cout<<"                              Default is to have as many colours as there are"<<std::endl;
    std::cout<<"                              different values for the coverage of Nodes in the graph"<<std::endl;
    std::cout<<"-c --graphColour     <TYPE>   Defines the range of colours to use for the output graph"<<std::endl;
    std::cout<<"                              the different types available are:"<<std::endl;
    std::cout<<"                              red-blue, blue-red, green-red-blue, red-blue-green"<<std::endl;
    std::cout<<"-L --longDescription          Set if you want the spacer sequence printed along with the ID in the spacer graph. [Default: false]"<<std::endl;
    std::cout<<"-G --showSingltons            Set if you want to print singleton spacers in the spacer graph [Default: false]"<<std::endl;
#ifdef DEBUG
    std::cout<<"--noDebugGraph                Stops creation of debug .gv files even if the DEBUG preprocessor macro is set [Default: false]"<<std::endl;
#endif
#ifdef RENDERING
    std::cout<<"--noRendering                 Stops rendering of .gv files even if the RENDERING preprocessor macro is set [Default: false]"<<std::endl;
#endif
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

//make a new directory
// Stoled from SaSSY
void recursiveMkdir(std::string dir) 
{
    std::string tmp;
    size_t pos = 0;
    while ( std::string::npos != (pos = dir.find('/',pos+1)) ) {
        tmp = dir.substr(0,pos);
        mkdir(tmp.c_str(), S_IRWXU);
    }
    mkdir(dir.c_str(), S_IRWXU);
}

int processOptions(int argc, char *argv[], options *opts) 
{
    int c;
    int index;
    bool scalling = false;
    while( (c = getopt_long(argc, argv, "a:b:c:d:D:f:Ghk:K:l:Ln:o:rs:S:Vw:x:y:", long_options, &index)) != -1 ) 
    {
        switch(c) 
        {
            case 'a':
#if RENDERING                
                if (strcmp(optarg, "dot") && strcmp(optarg, "neato") && strcmp(optarg, "fdp") && strcmp(optarg, "sfdp") && strcmp(optarg, "twopi") && strcmp(optarg, "circo")) {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: '"<<optarg<<"' is not a recognised layout algorithm. Please choose from the following:"<<std::endl;
#ifdef HAVE_NEATO
                    std::cerr<<"\tneato"<<std::endl;
#endif
#ifdef HAVE_DOT
                    std::cerr<<"\tdot"<<std::endl;
#endif
#ifdef HAVE_FDP
                    std::cerr<<"\tfdp"<<std::endl;
#endif
#ifdef HAVE_SFDP
                    std::cerr<<"\tsfdp"<<std::endl;
#endif
#ifdef HAVE_CIRCO
                    std::cerr<<"\tcirco"<<std::endl;
#endif
#ifdef HAVE_TWOPI
                    std::cerr<<"\ttwopi"<<std::endl;
#endif
                    
                } else {
                    opts->layoutAlgorithm = optarg;
                }
#else 
                std::cerr<<PACKAGE_NAME<<" [WARNING]: Not compiled with the correct options to allow rendering of graphs.\nEither use the --enable-rendering option during ./configure and make sure that the graphviz layout programs are in your PATH"<<std::endl;
#endif
                break;
            case 'b': 
                from_string<int>(opts->coverageBins, optarg, std::dec);
                if (opts->coverageBins < 1) 
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: The number of bins for coverage cannot be less than 1"<<std::endl;
                    usage();
                    exit(1);
                }
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
            case 'd': 
                from_string<unsigned int>(opts->lowDRsize, optarg, std::dec);
                if (opts->lowDRsize < 8) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for direct repeat sizes cannot be "<<opts->lowDRsize<<" changing to "<<CRASS_DEF_MIN_DR_SIZE<<std::endl;
                    opts->lowDRsize = CRASS_DEF_MIN_DR_SIZE;
                }
                break;
            case 'D': 
                from_string<unsigned int>(opts->highDRsize, optarg, std::dec);
                break;
            case 'f':
                from_string<int>(opts->covCutoff, optarg, std::dec);
                break;
            case 'G':
            	opts->showSingles= true;
                break;
            case 'h': 
                versionInfo(); 
                usage();
                exit(1); 
                break;
            case 'k': 
                from_string<int>(opts->kmer_clust_size, optarg, std::dec);
                if (opts->kmer_clust_size < 4) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Minimum value for kmer clustering size is: "<<4<<" changing to "<<CRASS_DEF_K_CLUST_MIN<<std::endl;
                    opts->kmer_clust_size = CRASS_DEF_K_CLUST_MIN;
                }
                break;
            case 'K': 
                from_string<int>(opts->cNodeKmerLength, optarg, std::dec);
                break;
            case 'l': 
                from_string<int>(opts->logLevel, optarg, std::dec);
                if(opts->logLevel > CRASS_DEF_MAX_LOGGING)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified log level higher than max. Changing log level to "<<CRASS_DEF_MAX_LOGGING<<" instead of "<<opts->logLevel<<std::endl;
                    opts->logLevel = CRASS_DEF_MAX_LOGGING;
                }
                break;
            case 'L':
                opts->longDescription = true;
                break;
            case 'n':
                from_string<unsigned int>(opts->minNumRepeats, optarg, std::dec);
                if (opts->minNumRepeats < 2) 
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: The mininum number of repeats cannot be less than 2"<<std::endl;
                    usage();
                    exit(1);
                }
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
                if (0 != stat(opts->output_fastq.c_str(),&file_stats)) 
                {
                    recursiveMkdir(opts->output_fastq);
                }
                // check that the directory is writable
                else if(access(optarg, W_OK))
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: You do not have permission to write to "<<optarg<<std::endl;
                    exit(1);
                }
                break;
            case 'r': 
                opts->reportStats = true; 
                break;
            case 's': 
                from_string<unsigned int>(opts->lowSpacerSize, optarg, std::dec);
                if (opts->lowSpacerSize < 8) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for spacer sizes cannot be "<<opts->lowSpacerSize<<" changing to "<<CRASS_DEF_MIN_SPACER_SIZE<<std::endl;
                    opts->lowSpacerSize = CRASS_DEF_MIN_SPACER_SIZE;
                }
                break;
            case 'S': 
                from_string<unsigned int>(opts->highSpacerSize, optarg, std::dec);
                break;
            case 'V': 
                versionInfo(); 
                exit(1); 
                break;
            case 'w': 
                from_string<unsigned int>(opts->searchWindowLength, optarg, std::dec);
                if ((opts->searchWindowLength < CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH) || (opts->searchWindowLength > CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH))
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified window length higher than max. Changing window length to " << CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH << " instead of " << opts->searchWindowLength<<std::endl;
                    // Change window length
                    opts->searchWindowLength = CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH;
                }
                break;        
            case 'x':
                from_string<double>(opts->averageSpacerScalling, optarg, std::dec);
                if (isNotDecimal(opts->averageSpacerScalling)) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The average spacer scalling must be a decimal number. Changing to "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<" instead of "<<opts->averageSpacerScalling<<std::endl;
                    opts->averageSpacerScalling = CRASS_DEF_HOMOPOLYMER_SCALLING;
                }
                scalling = true;
                break;
            case 'y':
                from_string<double>(opts->averageDrScalling, optarg, std::dec);
                if (isNotDecimal(opts->averageDrScalling)) 
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The average spacer scalling must be a decimal number. Changing to "<<CRASS_DEF_HOMOPOLYMER_SCALLING<<" instead of "<<opts->averageDrScalling<<std::endl;
                    opts->averageDrScalling = CRASS_DEF_HOMOPOLYMER_SCALLING;
                }
                scalling = true;
                break;
            case 0:
                if ( strcmp( "removeHomopolymers", long_options[index].name ) == 0 ) opts->removeHomopolymers = true;
                else if ( strcmp( "logToScreen", long_options[index].name ) == 0 ) opts->logToScreen = true;
                else if ( strcmp( "noScalling", long_options[index].name ) == 0 ) opts->dontPerformScalling = true;
#ifdef RENDERING 
                else if (strcmp( "noRendering", long_options[index].name ) == 0 ) opts->noRendering = true; 
#endif  
#ifdef DEBUG
                else if (strcmp( "noDebugGraph", long_options[index].name ) == 0 ) opts->noDebugGraph = true;
#endif
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
#ifdef PERFORM_CRASS_ASSEMBLY
   
    if (strcmp(argv[1], "assemble") == 0) 
    {
        // our user wants to do an assembly so let's load up the assembler main function
        return assemblyMain(argc - 1 , argv + 1);
    }
#endif
    /* application of default options */
    options opts;
    opts.logLevel              = CRASS_DEF_DEFAULT_LOGGING;              // level of verbosity allowed in the log file
    opts.reportStats           = CRASS_DEF_STATS_REPORT;                 // print a starts report currently not used
    opts.lowDRsize             = CRASS_DEF_MIN_DR_SIZE;                  // the lower size limit for a direct repeat
    opts.highDRsize            = CRASS_DEF_MAX_DR_SIZE;                  // the upper size limit for a direct repeat
    opts.lowSpacerSize         = CRASS_DEF_MIN_SPACER_SIZE;              // the lower limit for a spacer
    opts.highSpacerSize        = CRASS_DEF_MAX_SPACER_SIZE;              // the upper size limit for a spacer
    opts.output_fastq          = CRASS_DEF_OUTPUT_DIR;                   // the output directory for the output files
    opts.delim                 = CRASS_DEF_STATS_REPORT_DELIM;           // delimiter used in stats report currently not used
    opts.kmer_clust_size             = CRASS_DEF_K_CLUST_MIN;                  // number of kmers needed to be shared to add to a cluser
    opts.searchWindowLength    = CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH; // option 'w'used in long read search only
    opts.minNumRepeats         = CRASS_DEF_DEFAULT_MIN_NUM_REPEATS;      // option 'n'used in long read search only
    opts.logToScreen           = CRASS_DEF_LOGTOSCREEN;                  // log to std::cout rather than to the log file
    opts.removeHomopolymers    = CRASS_DEF_REMOVE_HOMOPOLYMERS;          // correct for homopolymer errors
    opts.coverageBins          = CRASS_DEF_NUM_OF_BINS;                  // The number of bins of colours
    opts.graphColourType       = CRASS_DEF_GRAPH_COLOUR;                 // the colour type of the graph
    opts.averageSpacerScalling = CRASS_DEF_HOMOPOLYMER_SCALLING;         // decimal for reduction in the spacer size
    opts.averageDrScalling     = CRASS_DEF_HOMOPOLYMER_SCALLING;         // decimal for the reduction in the direct repeat size
    opts.dontPerformScalling   = CRASS_DEF_NO_SCALLING;                  // turn all scalling off for the user to define variables
    opts.longDescription       = CRASS_DEF_SPACER_LONG_DESC;             // print a long description for the final spacer graph
    opts.showSingles           = CRASS_DEF_SPACER_SHOW_SINGLES;          // print singletons when making the spacer graph
    opts.cNodeKmerLength       = CRASS_DEF_NODE_KMER_SIZE;               // length of the kmers making up a crisprnode
    #ifdef DEBUG
    opts.noDebugGraph          = false;                                  // Even if DEBUG preprocessor macro is set do not produce debug graph files
    #endif
    #ifdef RENDERING
    opts.layoutAlgorithm       = DEFAULT_RENDERING_ALGORITHM;            // the graphviz layout algorithm to use
    opts.noRendering           = false;                                  // Even if RENDERING preprocessor macro is set do not produce any rendered images
    #else
    opts.layoutAlgorithm       = "unset";
    #endif
    opts.covCutoff             = CRASS_DEF_COVCUTOFF;

    int opt_idx = processOptions(argc, argv, &opts);

    if (opt_idx >= argc) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: Specify sequence files to process!"<<std::endl;
        usage();
        exit(1);
    }
    
    // always make a timestamp
    std::string logFileName;
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    
    strftime (buffer,80,"%d_%m_%Y_%H%M%S",timeinfo);
    std::string timestamp(buffer);

    // initialize the log file
    if (opts.logToScreen) 
    {
        logFileName = "";
    } 
    else 
    {
        // create a time stamp for the log file
        logFileName = opts.output_fastq + PACKAGE_NAME +"."+timestamp+".log";
        
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
    
    std::string cmd_line;
    for (int i = 0; i < argc; ++i) {
        cmd_line += argv[i];
        cmd_line += ' ';
    }
    
    WorkHorse * mHorse = new WorkHorse(&opts, timestamp,cmd_line);
    int error_code = mHorse->doWork(seq_files);
    delete mHorse;
    return error_code;
}
