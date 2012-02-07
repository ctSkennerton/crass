/*
 *  CrassException.h is part of the CRisprASSembler project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2012 Connor Skennerton & Michael Imelfort. All rights reserved. 
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

#ifndef crass_CrassException_h
#define crass_CrassException_h

#include <exception>
#include <sstream>
#include <string>
#include <cstring>

namespace crass {
    
    class exception 
    {
    public:
        // constructor
        exception(){}
        exception(const char * file, int line, const char * function ,const char * message)
        {
            std::stringstream ss;
            ss<<"[ERROR]: ";
            ss<< message<<std::endl;
            ss<<file<<" : "<<line<<" : "<<function;
            errorMsg = ss.str();
        }
        // destructor
        ~exception(){}
        
        virtual std::string what(void)
        {
            return errorMsg;
        }
    protected:
        std::string errorMsg;
        
    };
    
    class substring_exception: public exception{
    public:
        substring_exception(const char * message, const char * base_str, int start_pos, int length, const char * file, int line, const char * function)
        {
            std::stringstream ss;
            ss<<"[ERROR]: ";
            ss<< message<<std::endl;
            ss<<file<<" : "<<line<<" : "<<function<<std::endl;
            ss<<"Cutting substr: "<<start_pos<<" : "<<length<<std::endl;
            ss<<"From Read: "<<base_str<<std::endl;
            ss<<"Length: "<<strlen(base_str)<<std::endl;
            errorMsg = ss.str();

        }
        // destructor
        ~substring_exception(){}
        
        virtual std::string what(void)
        {
            return errorMsg;
        }
    protected:
        std::string errorMsg;
    };
}


#endif
