// File: SmithWaterman.h
// Original Author: Connor Skennerton
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Modified smithwaterman implementation
// 
// --------------------------------------------------------------------
//  Copyright  2011 Michael Imelfort and Connor Skennerton
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//
// --------------------------------------------------------------------

//////////////////////////////////////////////
// Simple Ends-Free Smith-Waterman Algorithm
//
// You will be prompted for input sequences
// Penalties and match scores are hard-coded
//
// Program does not perform multiple tracebacks if 
// it finds several alignments with the same score
//
// By Nikhil Gopal
// Similar implementation here: https://wiki.uni-koeln.de/biologicalphysics/index.php/Implementation_of_the_Smith-Waterman_local_alignment_algorithm
//////////////////////////////////////////////


#ifndef __SMITH_WATERMAN_H
  #define __SMITH_WATERMAN_H

// system includes
#include <string>
#include <map>

#define SW_MATCH                1
#define SW_MISMATCH             (-1)
#define SW_GAP                  (-1)
#define SW_SIM_SCORE(_a, _b)    ((_a == _b) ? SW_MATCH : SW_MISMATCH )

typedef std::pair<std::string, std::string> stringPair;

double findMax(double a, double b, double c, double d, int * index);

stringPair smithWaterman(std::string seqA, std::string seqB);
stringPair smithWaterman(std::string seqA, std::string seqB, int * aStartAlign, int * aEndAlign, int aStartSearch, int aEndSearch);

#endif // __SMITH_WATERMAN_H
