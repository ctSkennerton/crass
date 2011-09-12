// File: Rainbow.cpp
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Implementation of Rainbow functions
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
#include <math.h>
#include <stdio.h>

// local includes
#include "Rainbow.h"

// vertical scalers
#define RB_lower_offset (0.5)
#define RB_divisor (0.6666666666)
#define getValue(vAL) ((cos(vAL) + RB_lower_offset) * RB_divisor)

// crude rounding function

// setting initial limits
void Rainbow::setLimits(double lb, double ub, int res)
{
    //-----
    // set the limits
    //
    mLowerBound = lb;
    mUpperBound = ub;
    mResolution = res;
    
    // set the shortcut constants
    mScaleMultiplier = (mUpperScale - mLowerScale)/(mUpperBound - mLowerBound);
    mTickSize = (mUpperBound - mLowerBound)/(mResolution - 1);
}

void Rainbow::setType(RB_TYPE type)
{
    //-----
    // set the offsets based on type
    //
    mType = type;
    
    switch(type)
    {
        case RED_BLUE:
        {
            mRedOffset = 0;
            mGreenOffset = RB_divisor * PI * 2;
            mBlueOffset = RB_divisor * PI;
            
            mIgnoreRed = false;
            mIgnoreGreen = true;
            mIgnoreBlue = false;
            
            mLowerScale = 0;
            mUpperScale = (RB_divisor * PI);
            break;    
        }
        case RED_BLUE_GREEN:
        {
            mRedOffset = 0;
            mGreenOffset = RB_divisor * PI * 2;
            mBlueOffset = RB_divisor * PI;
            
            mIgnoreRed = false;
            mIgnoreGreen = false;
            mIgnoreBlue = false;
            
            mLowerScale = 0;
            mUpperScale = (RB_divisor * PI * 2);
            break;
        }
        case GREEN_BLUE_RED:
        {
            mRedOffset = RB_divisor * PI * 2;
            mGreenOffset = 0;
            mBlueOffset = RB_divisor * PI;

            mIgnoreRed = false;
            mIgnoreGreen = false;
            mIgnoreBlue = false;
            
            mLowerScale = 0;
            mUpperScale = (RB_divisor * PI * 2);
            break;    
        }
        case BLUE_RED:
        default:
        {
            // BLUE_RED,
            mRedOffset = RB_divisor * PI;
            mGreenOffset = RB_divisor * PI * 2;
            mBlueOffset = 0;
            
            mIgnoreRed = false;
            mIgnoreGreen = true;
            mIgnoreBlue = false;

            mLowerScale = 0;
            mUpperScale = (RB_divisor * PI);
        }
    }
    
    // reset this
    mScaleMultiplier = (mUpperScale - mLowerScale)/(mUpperBound - mLowerBound);
}

void Rainbow::setDefaults(void)
{
    //-----
    // Just waht it says
    //
    setType(RB_DEFAULT_TYPE);
    setLimits(RB_LB, RB_UB, RB_TICKS);
} 

// get a colour!;
std::string Rainbow::getColour(double value)
{
    //-----
    // Return a colour for the given value.
    // If nothing makes sense. return black
    //

    // check to see that everything makes sense    
    if(-1 == mResolution)
        return RB_ERROR_COLOUR;
    if(value > mUpperBound || value < mLowerBound)
        return RB_ERROR_COLOUR;
    
    // normalise the value to suit the ticks
    double normalised_value = round(value/mTickSize) * mTickSize;

    // map the normalised value onto the horizontal scale
    double scaled_value = (normalised_value - mLowerBound) * mScaleMultiplier + mLowerScale;
        
    // now construct the colour    
    std::stringstream ss;
    if(mIgnoreRed)    
        ss << "00";
    else
        ss << int2RGB(round(getValue(scaled_value - mRedOffset) * 255));

    if(mIgnoreGreen)    
        ss << "00";
    else
        ss << int2RGB(round(getValue(scaled_value - mGreenOffset) * 255));

    if(mIgnoreBlue)    
        ss << "00";
    else
        ss << int2RGB(round(getValue(scaled_value - mBlueOffset) * 255));

    return ss.str();
}

std::string Rainbow::int2RGB(int rgb)
{
    //----
    // convert an int between 0 and 255 to a hex RGB value
    //
    if(0 >= rgb)
    {
        return "00";
    }
    else
    {
        std::stringstream ss;
        if (16 > rgb)
            ss << "0" << std::hex << rgb;
        else
            ss << std::hex << rgb;
        return ss.str();
    }
}

