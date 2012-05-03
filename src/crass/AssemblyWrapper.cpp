/*
 *  AssemblyWrapper.cpp is part of the crass project
 *  
 *  Created by Connor Skennerton on 27/10/11.
 *  Copyright 2011 Connor Skennerton and Michael Imelfort. All rights reserved. 
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

// System Includes
#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

// Local Includes
#include "AssemblyWrapper.h"
#include "config.h"
#include "crassDefines.h"
#include "StlExt.h"
#include "Xml.h"
#include "kseq.h"
#include "SeqUtils.h"


KSEQ_INIT(gzFile, gzread)

void assemblyVersionInfo(void) 
{
    std::cout<<std::endl<<PACKAGE_FULL_NAME<<" ("<<PACKAGE_NAME<<")"<<std::endl<<"version "<<PACKAGE_MAJOR_VERSION<<" subversion "<<PACKAGE_MINOR_VERSION<<" revison "<<PACKAGE_REVISION<<" ("<<PACKAGE_VERSION<<")"<<std::endl<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Copyright (C) 2011 Connor Skennerton & Michael Imelfort"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY"<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it"<<std::endl;
    std::cout<<"under certain conditions: See the source for more details"<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
}

void assemblyUsage(void)
{
    std::cout<<"Usage: "PACKAGE_NAME<<" assemble ASSEMBLER -g INT -s LIST -x CRASS_XML_FILE -i INDIR [options]"<<std::endl;
    std::cout<<"\twhere ASSEMBLER is one of the assembly algorithms listed below:"<<std::endl<<std::endl;
#ifdef HAVE_VELVET
    std::cout<<"\tvelvet"<<std::endl;
#endif
    
#ifdef HAVE_CAP3
    std::cout<<"\tcap3"<<std::endl;
#endif
    std::cout<<std::endl;
    std::cout<< "-h --help                    This help message"<<std::endl;
    std::cout<< "-V --version                 Program and version information"<<std::endl;
#if 0
    std::cout<< "-l --logLevel        <INT>   Output a log file and set a log level [1 - "<<CRASS_DEF_MAX_LOGGING<<"]"<<std::endl;
    std::cout<< "--logToScreen                Print the logging information to screen rather than a file"<<std::endl;
#endif
    std::cout<< "-g --group           <INT>   ID of the group that you want to assemble"<<std::endl;
    std::cout<< "-s --segments        <LIST>  A comma separated list of numbered segments to assemble from the specified group"<<std::endl;
    std::cout<< "-x --xml             <FILE>  xml output file created by crass"<<std::endl;
    std::cout<< "-i --inDir           <DIR>   input directory for the assembly [default: .]"<<std::endl;    
#if 0
    std::cout<< "-p --pairedEnd               Set if you want paired end assembly.  The crass assembler will check for unmatched"<<std::endl;
    std::cout<<"                              pairs and add them into the input file."<<std::endl;
    std::cout<< "-I --insertSize      <INT>   size of the insert for paired end assembly"<<std::endl;
#endif
    std::cout<< "-o --outDir          <DIR>   name of the directory for the assembly output files"<<std::endl;



}

int processAssemblyOptions(int argc, char * argv[], assemblyOptions& opts)
{

	int c;
    int index;
    static struct option assemblyLongOptions [] = {
        
        {"segments", required_argument, NULL, 's'},
        {"group", required_argument, NULL, 'g'},
        {"logLevel", required_argument, NULL, 'l'},
        {"version", no_argument, NULL, 'V'},
        {"inDir", required_argument, NULL, 'i'},
        {"outDir", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'},
        {"pairedEnd", no_argument, NULL, 'p'},
        {"insertSize", required_argument, NULL, 'I'},
        {"logToScreen", no_argument, NULL, 0},
        {"xml",required_argument,NULL,'x'},
        {NULL, no_argument, NULL, 0}
    };
    try 
    {
        while( (c = getopt_long(argc, argv, "g:hi:I:l:o:ps:Vx:", assemblyLongOptions, &index)) != -1 ) 
        {
            switch(c) 
            {
                case 'g':
                {
                    from_string<int>(opts.group, optarg, std::dec);
                    break;
                }
                case 'h':
                {
                    assemblyUsage();
                    exit(1);
                    break;
                }
                case 'i':
                {
                    // Test to see if the file is ok.
                    struct stat inputDirStatus;
                    int iStat = stat(optarg, &inputDirStatus);
                    // stat failed
                    switch (iStat) 
                    {
                        case -1:
                        {
                            switch (errno)
                            {
                                case ENOENT:
                                {
                                    throw ( std::runtime_error("Input directory path does not exist, or path is an empty string.") );
                                    break;
                                }
                                case ELOOP:
                                {
                                    throw ( std::runtime_error("Too many symbolic links encountered while traversing the input directory path."));
                                    break;
                                }
                                case EACCES:
                                {
                                    throw ( std::runtime_error("You do not have permission to access the input directory."));
                                    break;
                                }
                                case ENOTDIR:
                                {
                                    throw ( std::runtime_error("Input is not a directory\n"));
                                    break;
                                }
                                default:
                                {
                                    throw (std::runtime_error("An error occured when reading the input directory"));
                                    break;
                                }
                            }
                            break;
                        }
                        default:
                        {
                            opts.inputDirName = optarg;
                            break;
                        }
                    }
                    break;
                }
                case 'I':
                {
                    from_string<int>(opts.insertSize, optarg, std::dec);
                    break;
                }
                case 'l':
                {
                    from_string<int>(opts.logLevel, optarg, std::dec);
                    if(opts.logLevel > CRASS_DEF_MAX_LOGGING)
                    {
                        std::cerr<<PACKAGE_NAME<<" [WARNING]: Specified log level higher than max. Changing log level to "<<CRASS_DEF_MAX_LOGGING<<" instead of "<<opts.logLevel<<std::endl;
                        opts.logLevel = CRASS_DEF_MAX_LOGGING;
                    }
                    break;
                }
                case 'o':
                {
                    // Test to see if the file is ok.
                    struct stat outputDirStatus;
                    int oStat = stat(optarg, &outputDirStatus);
                    switch (oStat) 
                    {
                        case -1:
                        {
                            switch (errno)
                            {
                                case ENOENT:
                                {
                                    RecursiveMkdir(optarg);
                                    break;
                                }
                                case ELOOP:
                                {
                                    throw ( std::runtime_error("Too many symbolic links encountered while traversing the output directory path."));
                                    break;
                                }
                                case EACCES:
                                {
                                    throw ( std::runtime_error("You do not have permission to access the output directory."));
                                    break;
                                }
                                case ENOTDIR:
                                {
                                    throw ( std::runtime_error("Output path is not a directory\n"));
                                    break;
                                }
                                default:
                                {
                                    throw std::runtime_error("An error occurred when accessing the output directory");
                                    break;
                                }
                            }
                        }
                        default:
                        {
                            opts.outputDirName = optarg;
                            break;
                        }
                    }
                    break;
                }
                case 'p':
                {
                    opts.pairedEnd = true;
                    break;
                }
                case 's':
                {
                    opts.segments = optarg;
                    break;
                }
                case 'V':
                {
                    assemblyVersionInfo();
                    exit(1);
                    break;
                }
                case 'x':
                {
                    // Test to see if the file is ok.
                    struct stat inputDirStatus;
                    int xStat = stat(optarg, &inputDirStatus);
                    // stat failed
                    switch (xStat) 
                    {
                        case -1:
                        {
                            switch (errno)
                            {
                                case ENOENT:
                                {
                                    throw ( std::runtime_error("Path to the XML file does not exist, or path is an empty string.") );
                                    break;
                                }
                                case ELOOP:
                                {
                                    throw ( std::runtime_error("Too many symbolic links encountered while traversing the path to the XML file."));
                                    break;
                                }
                                case EACCES:
                                {
                                    throw ( std::runtime_error("You do not have permission to access the XML file."));
                                    break;
                                }
                                default:
                                {
                                    throw (std::runtime_error("An error occured when reading the XML file"));
                                    break;
                                }
                            }
                            break;
                        }
                        default:
                        {
                            opts.xmlFileName = optarg;
                            break;
                        }
                    }
                    break;
                }
                case 0:
                {
                    if ( strcmp( "logToScreen", assemblyLongOptions[index].name ) == 0 ) opts.logToScreen = true;
                    break;
                }
                default:
                {
                    assemblyUsage();
                    exit(1);
                    break;
                }
            }
        }
        // make sure that the required flags are set
        if (opts.xmlFileName.empty()) 
        {
            throw (std::runtime_error("You must specify an xml file with -x"));
        }
        if (!(opts.group)) 
        {
            throw (std::runtime_error("You must specify a group number with -g"));
        }
        if (opts.segments.empty()) 
        {
            throw (std::runtime_error("You must specify a list of contigs with -s"));
        }
        if (opts.inputDirName.empty()) 
        {
            throw (std::runtime_error("You must specify an input directory with -i"));
        }

    } 
    catch (std::exception& e) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: "<<e.what()<<std::endl;
        exit(1);
    }

    
    return optind;
}


void parseSegmentString(std::string& segmentString, std::set<std::string>& segments)
{
    // split the segment id string into the individual pieces
    std::vector<std::string> tmp;
    tokenize(segmentString, tmp, ",");

    std::vector<std::string>::iterator tmp_iter = tmp.begin();
    while (tmp_iter != tmp.end()) 
    {
        // prepend on a 'C' cause that is the format of the xml file
        segments.insert("C" + *tmp_iter);

        tmp_iter++;
    }
}

void generateTmpAssemblyFile(std::string fileName, std::set<std::string>& wantedContigs, assemblyOptions& opts, std::string& tmpFileName)
{
    
    gzFile fp = getFileHandle((opts.inputDirName + fileName).c_str());
    kseq_t *seq;
    int l;
    
    // initialize an output file handle
    std::ofstream out_file;
    out_file.open(tmpFileName.c_str());
    
    if (out_file.good()) 
    {
        // initialize seq
        seq = kseq_init(fp);
        
        // read sequence  
        while ( (l = kseq_read(seq)) >= 0 ) 
        {
            std::string name = seq->name.s;
            // get the last char of the name and see if it's from one of our contigs
            std::stringstream ss;
            ss << name.substr(name.length() - 2);
            std::string contig_number = ss.str();
            if (wantedContigs.find(contig_number) != wantedContigs.end()) 
            {
                // this read comes from a segment that we want
                
                // check to see if it is fasta or fastq
                if (seq->qual.s) 
                {
                    // it's fastq
                    out_file<<'@'<<seq->name.s<<std::endl;
                    out_file<<seq->seq.s<<std::endl;
                    out_file<<'+';
                    if(seq->comment.s) 
                    {
                        out_file<<seq->comment.s;
                    }
                    out_file<<std::endl<<seq->qual.s<<std::endl;
                } 
                else 
                {
                    // it's fasta
                    out_file<<'>'<<seq->name.s;
                    if (seq->comment.s) 
                    {
                        out_file<<' '<<seq->comment.s;
                    }
                    out_file<<std::endl<<seq->seq.s<<std::endl;
                }                
            }
        }
    }
}


int velvetWrapper( int hashLength, assemblyOptions& opts, std::string& tmpFileName)
{
    // get the hash length from the DR
    
    try 
    {
        // create the command string for velvet
        std::string h_cmd = "velveth " + opts.outputDirName +" " + to_string(hashLength) + " " + opts.inputDirName + tmpFileName;
        std::string g_cmd = "velvetg " + opts.outputDirName;     
        int h_exit = system(h_cmd.c_str());
        if (h_exit) 
        {
            // print error
            throw (std::runtime_error("velveth did not exit normally"));
        }
        int g_exit = system(g_cmd.c_str());
        if (g_exit) 
        {
            // print error
            throw (std::runtime_error("velvetg did not exit normally"));
        }
    } 
    catch (std::exception& e) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: "<<e.what()<<std::endl;
        return 1;
    }  
    return 0;
    
}

int capWrapper(int overlapLength, assemblyOptions& opts, std::string& tmpFileName)
{
    // cap3 doesn't control the output directory so we need to move 
    // the tmp_file to the output directory
    try 
    {
        std::ifstream input_file;
        input_file.open((opts.inputDirName + tmpFileName).c_str());
        if (input_file.good()) 
        {
            std::ofstream output_file;
            output_file.open((opts.outputDirName + tmpFileName).c_str());
            if (output_file.good()) 
            {
                char c;
                while (input_file.get(c)) output_file << c;
            } 
            else 
            {
                throw (std::runtime_error("Output file stream is not good"));
            }
        } 
        else 
        {
            throw (std::runtime_error("Input file stream is not good"));
        }
        std::string cap3cmd = "cap3 " + opts.outputDirName + tmpFileName + " -o " + to_string(overlapLength) + " -x crass > " + opts.outputDirName + tmpFileName + ".crass.cap3";
        int cap_exit = system(cap3cmd.c_str());
        if (cap_exit) 
        {
            throw (std::runtime_error("cap3 did not exit normally"));
        }
    } 
    catch (std::exception& e) 
    {
        std::cerr<<PACKAGE_NAME<<" [ERROR]: "<<e.what()<<std::endl;
        return 1;
    }
    return 0;
}


int assemblyMain(int argc, char * argv[])
{
    if (argc == 1) 
    {
        // print assembly wrapper help
        assemblyUsage();
        exit(1);
    }
    
    assemblyOptions opts;
    ASSEMBLERS wantedAssembler = unset;

#ifdef HAVE_VELVET
    if (strcmp(argv[1], "velvet") == 0) 
    {
        // the user wants velvet
        wantedAssembler = velvet;
    } 
#endif
    
#ifdef HAVE_CAP3
    if(strcmp(argv[1], "cap3") == 0)
    {
        // the user wants CAP3
        wantedAssembler = cap3;
    }
#endif
    
    if(unset == wantedAssembler)
    {
        std::cout << "**ERROR: No valid assemblers installed" << std::endl;
        return 43;
    }

    processAssemblyOptions(argc - 1, argv + 1, opts);
    
    std::set<std::string> segments;
    parseSegmentString(opts.segments, segments);
    
    std::list<std::string> spacers_for_assembly;

    //parse xml file
    crispr::xml::reader xml_parser;  
    std::string direct_repeat;
    std::string group_as_string = "G" + to_string(opts.group);

    xml_parser.parseXMLFile(opts.xmlFileName, group_as_string, &direct_repeat, segments, spacers_for_assembly );
    
    // parse the read file and create tmp.fa for the right segments
    //build the read file name from what we know
    std::string group_read_file = "Group_" + to_string(opts.group) + "_" + direct_repeat + ".fa";
    
    // get the tmp file name
    std::string tmp_file_name = PACKAGE_NAME;
    tmp_file_name += "_tmp.fa";
    
    generateTmpAssemblyFile(group_read_file, segments, opts, tmp_file_name);
    int return_value = 42;
    switch (wantedAssembler) 
    {
        case velvet:
            // velvet wrapper
            return_value = velvetWrapper(calculateOverlapLength((int)direct_repeat.length()), opts, tmp_file_name);
            break;
         case cap3:
            // cap3 wrapper
           return_value = capWrapper(calculateOverlapLength((int)direct_repeat.length()), opts, tmp_file_name);
            break;
        default:
            // assembler not known throw error
            return_value = 1;
            break;
    }
    return return_value;
}

