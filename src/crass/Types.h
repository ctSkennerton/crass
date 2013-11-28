/*
 *  Types.h is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011, 2012 Connor Skennerton & Michael Imelfort. All rights reserved. 
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

#ifndef crass_Types_h
#define crass_Types_h
#include <map>
#include <vector>
#include <string>
#include "ReadHolder.h"
#include "StringCheck.h"


// forward declaration of readholder class
//class ReadHolder;
namespace crass {
    typedef std::vector<std::string> list;
}

// Types cut from libcrispr.h
typedef std::map<std::string, bool> lookupTable;

typedef std::vector<ReadHolder *> ReadList;
typedef std::vector<ReadHolder *>::iterator ReadListIterator;

// direct repeat as a string and a list of the read objects that contain that direct repeat
typedef std::map<StringToken, ReadList *> ReadMap;
typedef std::map<StringToken, ReadList *>::iterator ReadMapIterator;

// Types from WorkHorse.h
// for storing clusters of DRs
// indexed using StringCheck type tokens
typedef std::vector<StringToken> DR_Cluster; 
typedef std::vector<StringToken>::iterator DR_ClusterIterator;

typedef std::map<int, DR_Cluster *>::iterator DR_Cluster_MapIterator;
typedef std::map<int, DR_Cluster *> DR_Cluster_Map;

typedef std::map<int, std::map<std::string, int> * > GroupKmerMap;

typedef std::vector<std::string> Vecstr;

#endif
