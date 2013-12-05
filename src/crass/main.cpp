/*
 *  main.cpp is part of the CRisprASSembler project
 *
 *  Created by Connor Skennerton.
 *  Copyright 2011 - 2013 Connor Skennerton & Michael Imelfort. All rights reserved.
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
#include "main.h"
#include "config.h"
#include "StlExt.h"
#include "Storage.h"
#include "Search.h"
#include <iostream>
#include <sys/stat.h>
#include <getopt.h>
#include <string>
#define OUT_FILE_PERMISSION_MASK (S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)

void RecursiveMkdir(std::string dir)
{
    std::string tmp;
    size_t pos = 0;
    while ( std::string::npos != (pos = dir.find('/',pos+1)) ) {
        tmp = dir.substr(0,pos);
        mkdir(tmp.c_str(), OUT_FILE_PERMISSION_MASK);
    }
    mkdir(dir.c_str(), OUT_FILE_PERMISSION_MASK);
}

void versionInfo(void)
{
    std::cout<<PACKAGE_NAME<<" version "<<PACKAGE_VERSION<<std::endl<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Copyright (C) 2011-2013 Connor Skennerton & Michael Imelfort"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY"<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it"<<std::endl;
    std::cout<<"under certain conditions: See the source for more details"<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
}

void usage() {
    std::cout<<std::endl;
    std::cout<<"Usage:  "<<std::endl;
    std::cout<<PACKAGE_NAME<<"  [options] <inputFile>..."<<std::endl;
    //std::cout<<PACKAGE_NAME<<"  [options] -ps (<readfile1> <readfile2>)..."<<std::endl;
    //std::cout<<PACKAGE_NAME<<"  [options] -pi <readfile12>..."<<std::endl;

    std::cout<<"Options:"<<std::endl;
    std::cout<< "-h --help                    This help message"<<std::endl;
    std::cout<< "-l --logLevel        <INT>   Output a log file and set a log level [1 - "<<CRASS_DEF_MAX_LOGGING<<"]"<<std::endl;
    std::cout<< "-L --logFile         <FILE>  Name of the file for log in the output directory."<<std::endl;
    std::cout<< "                             use 'stdout' or 'stderr' to print logging info to that stream"<<std::endl;
    std::cout<< "                             [Default: log.txt]"<<std::endl;
    std::cout<< "-o --outDir          <DIR>   Output directory; will be created if it doesn't exist [default: .]"<<std::endl;
    std::cout<< "-V --version                 Program and version information"<<std::endl;
    std::cout<< "-d --minDR           <INT>   Minimim length of the direct repeat"<<std::endl;
    std::cout<< "                             to search for [Default: "<<CRASS_DEF_MIN_DR_SIZE<<"]"<<std::endl;
    std::cout<< "-D --maxDR           <INT>   Maximim length of the direct repeat"<<std::endl;
    std::cout<< "                             to search for [Default: "<<CRASS_DEF_MAX_DR_SIZE<<"]"<<std::endl;
    std::cout<< "-s --minSpacer       <INT>   Minimim length of the spacer to search for [Default: "<<CRASS_DEF_MIN_SPACER_SIZE<<"]"<<std::endl;
    std::cout<< "-S --maxSpacer       <INT>   Maximim length of the spacer to search for [Default: "<<CRASS_DEF_MAX_SPACER_SIZE<<"]"<<std::endl;
    std::cout<< "-f --covCutoff       <INT>   Remove groups with less than x spacers [Default: "<<CRASS_DEF_COVCUTOFF<<"]"<<std::endl;
    std::cout<< "-k --kmerCount       <INT>   The number of the kmers that need to be"<<std::endl;
    std::cout<< "                             shared for clustering [Default: "<<CRASS_DEF_K_CLUST_MIN<<"]"<<std::endl;
    std::cout<<std::endl;


}

int processOptions(int argc, char *argv[], CrassOpts& opts)
{
    static struct option long_options [] = {
        
        {"minDR", required_argument, NULL, 'd'},
        {"maxDR", required_argument, NULL, 'D'},
        {"covCutoff",required_argument,NULL,'f'},
        {"help", no_argument, NULL, 'h'},
        {"kmerCount", required_argument, NULL, 'k'},
        {"logLevel", required_argument, NULL, 'l'},
        {"logFile", required_argument, NULL, 'L'},
        {"outDir", required_argument, NULL, 'o'},
        {"minSpacer", required_argument, NULL, 's'},
        {"maxSpacer", required_argument, NULL, 'S'},
        {"version", no_argument, NULL, 'V'},

        {NULL, no_argument, NULL, 0}
    };
    
    int c;
    int index;
    bool scalling = false;
    while( (c = getopt_long(argc, argv, "d:D:f:hk:l:L:o:s:S:V", long_options, &index)) != -1 )
    {
        switch(c)
        {
            case 'd':
                from_string<unsigned int>(opts.lowDRsize, optarg, std::dec);
                if (opts.lowDRsize < 8)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for direct repeat sizes cannot be "<<opts.lowDRsize<<" changing to "<<23<<std::endl;
                    opts.lowDRsize = 23;
                }
                break;
            case 'D':
                from_string<unsigned int>(opts.highDRsize, optarg, std::dec);
                break;
            case 'f':
                from_string<int>(opts.covCutoff, optarg, std::dec);
                break;
            case 'h':
                versionInfo();
                usage();
                exit(1);
                break;
            case 'k':
                from_string<int>(opts.kmerClustSize, optarg, std::dec);
                if (opts.kmerClustSize < 4)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Minimum value for kmer clustering size is: "<<4<<" changing to "<<6<<std::endl;
                    opts.kmerClustSize = 6;
                }
                break;
            case 'l':
                from_string<int>(opts.logLevel, optarg, std::dec);
                if(opts.logLevel > 10)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified log level higher than max. Changing log level to "<<10<<" instead of "<<opts.logLevel<<std::endl;
                    opts.logLevel = 10;
                }
                break;
            case 'L':
                opts.logFile = optarg;
            case 'o':
                opts.outdir = optarg;
                opts.outdir += '/';
                
                // check if our output folder exists
                struct stat file_stats;
                if (0 != stat(opts.outdir.c_str(),&file_stats))
                {
                    RecursiveMkdir(opts.outdir);
                }
                // check that the directory is writable
                else if(access(optarg, W_OK))
                {
                    std::cerr<<PACKAGE_NAME<<" [ERROR]: You do not have permission to write to "<<optarg<<std::endl;
                    exit(1);
                }
                break;
            case 's':
                from_string<unsigned int>(opts.lowSpacerSize, optarg, std::dec);
                if (opts.lowSpacerSize < 8)
                {
                    std::cerr<<PACKAGE_NAME<<" [WARNING]: The lower bound for spacer sizes cannot be "<<opts.lowSpacerSize<<" changing to "<<CRASS_DEF_MIN_SPACER_SIZE<<std::endl;
                    opts.lowSpacerSize = CRASS_DEF_MIN_SPACER_SIZE;
                }
                break;
            case 'S':
                from_string<unsigned int>(opts.highSpacerSize, optarg, std::dec);
                break;
            case 'V':
                versionInfo();
                exit(1);
                break;
            default:
                versionInfo();
                usage();
                exit(1);
                break;
        }
    }
    // Sanity checks for the high and low dr size
    if (opts.lowDRsize >= opts.highDRsize)
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: The lower direct repeat bound is bigger than the higher bound ("<<opts.lowDRsize<<" >= "<<opts.highDRsize<<")"<<std::endl;
        usage();
        exit(1);
    }
    // Sanity checks for the high and low spacer size
    if (opts.lowSpacerSize >= opts.highSpacerSize)
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: The lower spacer bound is bigger than the higher bound ("<<opts.lowSpacerSize<<" >= "<<opts.highSpacerSize<<")"<<std::endl;
        usage();
        exit(1);
    }
    
    return optind;
}


int main(int argc, char * argv[]) {
    if(argc == 1)
    {
        usage();
        return 1;
    }
    CrassOpts opts;
    int opt_idx = processOptions(argc, argv, opts);
    if(opts.logFile == "stderr" || opts.logFile == "stdout") {
        intialiseGlobalLogger(opts.logFile, opts.logLevel);
    } else {
        intialiseGlobalLogger(opts.outdir + opts.logFile, opts.logLevel);
    }
    logTimeStamp();
    if (opt_idx >= argc)
    {
        logError("no sequence files to process");
        std::cerr<<PACKAGE_NAME<<" [ERROR]: Specify sequence files to process!"<<std::endl;
        usage();
        return 1;
    }
    
    crass::Storage read_store = crass::Storage();
    crass::Search searcher = crass::Search(&read_store);
    for (; opt_idx < argc; ++opt_idx) {
        searcher.searchFileSerial(argv[opt_idx]);
    }
    read_store.clusterRepeats(opts.outdir, 0.9, 1);
    
    searcher.searchFileForPatterns(argv[opt_idx]);
    
    read_store.inspect(std::cout);
    return 0;
}