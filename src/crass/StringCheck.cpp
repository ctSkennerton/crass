// File: StringCheck.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Implementation of StringCheck functions
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
// system includes
#include <iostream>
#include <sstream>

// local includes
#include "StringCheck.h"
#include "LoggerSimp.h"


StringToken StringCheck::addString(std::string& newStr)
{
    //-----
    // add the string and retuen it's token
    //
    mNextFreeToken++;
    mT2S_map[mNextFreeToken] = newStr;
    mS2T_map[newStr] = mNextFreeToken;
    return mNextFreeToken;
}

StringToken StringCheck::addString(const char * newStr)
{
    //-----
    // add the string and retuen it's token
    //
    mNextFreeToken++;
    mT2S_map[mNextFreeToken] = newStr;
    mS2T_map[newStr] = mNextFreeToken;
    return mNextFreeToken;
}

std::string StringCheck::getString(const StringToken token)
{
    //-----
    // return the string for a given token or spew
    //
    if(mT2S_map.find(token) != mT2S_map.end()) {
        return mT2S_map[token];
    } else {
        logError("Token: "<<token<<" not stored!");
    }
    return "";
}

StringToken StringCheck::getToken( std::string& queryStr)
{
    //-----
    // return the token or 0
    //
    if(mS2T_map.end() == mS2T_map.find(queryStr))
        return 0;
    return mS2T_map[queryStr];
}

StringToken StringCheck::getToken( const char * queryStr)
{
    //-----
    // return the token or 0
    //
    if(mS2T_map.end() == mS2T_map.find(queryStr))
        return 0;
    return mS2T_map[queryStr];
}
