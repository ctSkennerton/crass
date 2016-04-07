/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * main.cpp
 * Copyright (C) Connor Skennerton 2011, 2012 <c.skennerton@gmail.com>
 * 
 * crisprtools is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * crisprtools is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <getopt.h>
#include <cstring>

#include "config.h"
#include "Utils.h"
#include "MergeTool.h"
#include "SplitTool.h"
#include "ExtractTool.h"
#include "FilterTool.h"
#include "SanitiseTool.h"
#if RENDERING && HAVE_LIBCDT && HAVE_LIBGRAPH && HAVE_LIBGVC
#include "DrawTool.h"
#endif
#include "StatTool.h"
#include "RemoveTool.h"
void usage (void)
{
	std::cout<<PACKAGE_NAME<<" ("<<PACKAGE_VERSION<<")"<<std::endl;
	std::cout<<PACKAGE_NAME<<" is a set of smal utilities for manipulating .crispr files"<<std::endl;
	std::cout<<"The .crispr file specification is a standard xml based format for describing CRISPRs"<<std::endl;
	std::cout<<"Type "<<PACKAGE_NAME<<" <subcommand> -h for help on each utility"<<std::endl;
	std::cout<<"Usage:\t"<<PACKAGE_NAME<<" <subcommand> [options]"<<std::endl<<std::endl;
    std::cout<<"subcommand:  merge       combine multiple files"<<std::endl;
	std::cout<<"             help        display this message and exit"<<std::endl;
	std::cout<<"             extract     extract sequences in fasta"<<std::endl;
	std::cout<<"             filter      make new files based on parameters"<<std::endl;
	std::cout<<"             sanitise    change the IDs of elements"<<std::endl;
#if RENDERING && HAVE_LIBCDT && HAVE_LIBGRAPH && HAVE_LIBGVC
    std::cout<<"             draw        create a rendered image of the CRISPR with Graphviz"<<std::endl;
#endif
	std::cout<<"             stat        show statistics on some or all CRISPRs"<<std::endl;
    std::cout<<"             rm          remove a group from a .crispr file"<<std::endl;
}

int main(int argc, char ** argv)
{
	if(argc == 1)
	{
		usage();
		return 1;
	}
	else if(!strcmp(argv[1], "help"))  {usage() ; return 0;}
	else if(!strcmp(argv[1], "merge")) return mergeMain(argc -1 , argv + 1);
    else if(!strcmp(argv[1], "split")) return splitMain(argc - 1, argv + 1);
	else if(!strcmp(argv[1], "extract")) return extractMain(argc - 1 , argv + 1);
	else if(!strcmp(argv[1],"filter")) return filterMain(argc - 1, argv + 1);
	else if(!strcmp(argv[1], "sanitise")) return sanitiseMain(argc - 1 , argv + 1);
#if RENDERING && HAVE_LIBCDT && HAVE_LIBGRAPH && HAVE_LIBGVC
	else if(!strcmp(argv[1], "draw")) return drawMain(argc - 1, argv + 1);
#endif
	else if(!strcmp(argv[1], "stat")) return statMain(argc - 1, argv + 1);
	else if (!strcmp(argv[1], "rm")) return removeMain(argc -1 , argv + 1);
	else
	{
		std::cerr<<"Unknown option: "<<argv[1]<<std::endl;
		usage();
	}
	return 0;
}
