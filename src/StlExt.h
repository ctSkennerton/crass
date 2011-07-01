// File: StlExt.h
// Original Author: Connor Skennerton on 1/07/11
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Place to put all ya'll stl extensions
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


template <class T1, class T2>
void addOrIncrement(std::map<T1, T2> &inMap, T1 &searchThing)
{
    
    if (inMap.find(searchThing) != inMap.end())
    {
        inMap[searchThing] += 1;
    }
    else
    {
        inMap[searchThing] = 1;
    }
}