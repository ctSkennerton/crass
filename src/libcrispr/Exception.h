/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * crisprtools
 * Copyright (C) Connor Skennerton 2011 <c.skennerton@gmail.com>
 * 
crisprtools is free software: you can redistribute it and/or modify it
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

#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <exception>
#include <string>
#include <sstream>
namespace crispr {
    
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
    
    class input_exception: public exception{
        public:
        input_exception(const char * message)
        {
            errorMsg = message;
        }
        // destructor
        ~input_exception(){}
        
        virtual std::string what(void)
        {
            return errorMsg;
        }
    protected:
        std::string errorMsg;
    };
    
    class xml_exception: public exception{
    public:
        xml_exception(const char * file, int line, const char * function ,const char * message)
        {
            std::stringstream ss;
            ss<<"[XML_ERROR]: ";
            ss<< message<<std::endl;
            ss<<file<<" : "<<line<<" : "<<function;
            errorMsg = ss.str();
        }
        // destructor
        ~xml_exception(){}
        
        virtual std::string what(void)
        {
            return errorMsg;
        }
    protected:
        std::string errorMsg;
    };
    
    class runtime_exception: public exception{
    public:
        runtime_exception(const char * file, int line, const char * function ,const char * message)
        {
            std::stringstream ss;
            ss<<"[RUNTIME_ERROR]: ";
            ss<< message<<std::endl;
            ss<<file<<" : "<<line<<" : "<<function;
            errorMsg = ss.str();
        }
        // destructor
        ~runtime_exception(){}
        
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
            ss<<"[SUBSTRING_ERROR]: ";
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
#endif // _CRISPREXCEPTION_H_
