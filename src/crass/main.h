/*
 *  main.h is part of the CRisprASSembler project
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

#ifndef __crass__main__
#define __crass__main__

#include "crassDefines.h"

struct CrassOpts {
    int                 logLevel;                                           // level of verbosity allowed in the log file
    unsigned int        lowDRsize;                                          // the lower size limit for a direct repeat
    unsigned int        highDRsize;                                         // the upper size limit for a direct repeat
    unsigned int        lowSpacerSize;                                      // the lower limit for a spacer
    unsigned int        highSpacerSize;                                     // the upper size limit for a spacer
    std::string         outdir;                                             // the output directory for the output files
    int                 kmerClustSize;                                    // number of kmers needed to be shared to add to a cluser
    std::string         logFile;                                            // Log name default = log.txt, use 'stdout' for screen
    int                 covCutoff;                                          // The lower bounds of acceptable numbers of reads that a group can have
    
    CrassOpts() :
        logLevel(1),
        lowDRsize(23),
        highDRsize(47),
        lowSpacerSize(26),
        highSpacerSize(50),
        outdir("./"),
        kmerClustSize(6),
        logFile("log.txt"),
        covCutoff(3)
    {}
    
};
#endif /* defined(__crass__main__) */
