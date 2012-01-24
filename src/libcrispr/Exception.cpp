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

#include "Exception.h"
#include <sstream>
namespace crispr {
    exception::exception (const char * file, int line, const char * function ,const char * message)
    {
        std::stringstream ss;
        ss<<"[ERROR]: ";
        ss<< message<<std::endl;
        ss<<file<<" : "<<line<<" : "<<function;
        errorMsg = ss.str();
    }

    input_exception::input_exception( const char * message)
    {
        errorMsg = message;
    }
    
    xml_exception::xml_exception( const char * file, int line, const char * function ,const char * message)
    {
        std::stringstream ss;
        ss<<"[XML_ERROR]: ";
        ss<< message<<std::endl;
        ss<<file<<" : "<<line<<" : "<<function;
        errorMsg = ss.str();
    }
    
    runtime_exception::runtime_exception( const char * file, int line, const char * function ,const char * message)
    {
        std::stringstream ss;
        ss<<"[RUNTIME_ERROR]: ";
        ss<< message<<std::endl;
        ss<<file<<" : "<<line<<" : "<<function;
        errorMsg = ss.str();
    }
}
