/*
 *  crass.h is part of the CRisprASSembler project
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

#ifndef __CRASS_H
    #define __CRASS_H

// system includes
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

// local includes
#include "crassDefines.h"
#include "kseq.h"
#include "CrisprNode.h"
#include "WorkHorse.h"

//**************************************
// user input + system
//**************************************

static struct option long_options [] = {

    {"minDR", required_argument, NULL, 'd'},
    {"maxDR", required_argument, NULL, 'D'},
    {"minSpacer", required_argument, NULL, 's'},
    {"maxSpacer", required_argument, NULL, 'S'},
    {"logLevel", required_argument, NULL, 'l'},
    {"maxMismatches", required_argument, NULL, 'm'},
    {"version", no_argument, NULL, 'V'},
    {"kmerCount", required_argument, NULL, 'k'},
    {"outDir", required_argument, NULL, 'o'},
    {"help", no_argument, NULL, 'h'},
    {"reportStats", no_argument, NULL, 'r'},
    {"dumpReads", no_argument, NULL, 0},
    {"windowLength", required_argument, NULL, 'w'},
    {"minNumRepeats", required_argument, NULL, 'n'},
    {"logToScreen", no_argument, NULL, 0},
    {"removeHomopolymers",no_argument,NULL,0},
    {NULL, no_argument, NULL, 0}
};

void  usage(void);

void  help(void);

void  versionInfo(void);

int   processOptions(int argc, char *argv[], options *opts);





#endif
