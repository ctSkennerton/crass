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
#include "crass_defines.h"
#include "kseq.h"
#include "CrisprNode.h"
#include "WorkHorse.h"

//**************************************
// user input + system
//**************************************

void  readsHelp(void);

void  mainHelp(void);

void  version_info(void);

int   processReadsOptions(int argc, char *argv[], options *opts);

int   genomeMain(int argc, char * argv[]);

int   readsMain(int argc, char * argv[]);

int   processGenomeOptions(int argc, char * argv[], options * opts);


#endif
