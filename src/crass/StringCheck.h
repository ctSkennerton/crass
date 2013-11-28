// File: StringCheck.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Give this guy a string, get a token, give this guy a token, get a string
// All token are unique, all strings aren't!
// 
// Basically a glorified map
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

#ifndef StringCheck_h
#define StringCheck_h

// system includes
#include <iostream>
#include <Unordered_map>


// typedefs
typedef int StringToken;

class StringCheck 
{
public:
    StringCheck(std::string name) : mNextFreeToken(1), mName(name){}
		StringCheck(void) : mNextFreeToken(1), mName("unset"){}
        ~StringCheck(void) {}  
        
        StringToken addString(std::string newStr);
        std::string getString(StringToken token);
        StringToken getToken(std::string queryStr);
        
        inline void setName(std::string name) { mName = name; }
private:
        // members
        StringToken mNextFreeToken;                            // der
        std::unordered_map<StringToken, std::string> mT2S_map;           // token to string map
        std::unordered_map<std::string, StringToken> mS2T_map;           // string to token map
        
        std::string mName;
};

#endif //StringCheck_h
