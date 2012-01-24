// File: Rainbow.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Simple heat map generator. Returns hex-RGB colours based on limits set by user
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

#ifndef Rainbow_h
#define Rainbow_h

// system includes
#include <iostream>
#include <math.h>

// local includes

// defaults
#define RB_ERROR_COLOUR "000000"    // error colour. return this whenever you want to panic
#define RB_DEFAULT_TYPE BLUE_RED    // default type of heatmap
#define RB_LB 0                     // default upper bound
#define RB_UB 1                     // default lowe bound
#define RB_TICKS 10                 // default ticks

// some math constants we need (I know!)
#define PI (3.1415927)

// we can have anumber of different types of heatmaps
enum RB_TYPE
{
    RED_BLUE,
    BLUE_RED,
    RED_BLUE_GREEN,
    GREEN_BLUE_RED
};

class Rainbow 
{
    public:
        Rainbow() { setDefaults(); } 
        ~Rainbow() {}  
        void setDefaults(void);
        
        // set the limits we need!
        void setLimits(double lb, double ub, int res);
        void setLimits(double lb, double ub) { setLimits(lb, ub, (int)(ub-lb)+1); }
        void setType(RB_TYPE type);
    
        void setUpperBound(double ub){setLimits(mLowerBound, ub, (int)(ub - mLowerBound)+1); }
        void setLowerBound(double lb){setLimits(lb, mUpperBound, (int)(mUpperBound - lb)+1); }

        // get a colour!
        std::string getColour(double value);
        std::string int2RGB(int rgb);
        
        // get the limits
        double getUpperLimit(void) { return mUpperBound; }
        double getLowerLimit(void) { return mLowerBound; }
        
    private:
        // members
        double mLowerBound;        // lowest number we'll colour`
        double mUpperBound;        // highest
        int mResolution;           // number of steps between low and high (inclusive)
        
        RB_TYPE mType;          // type of the heatmap 
        
        // We need offsets for the start of each cos graph
        double mRedOffset;
        double mGreenOffset;
        double mBlueOffset;
        
        // should we ignore any colour
        bool mIgnoreRed;
        bool mIgnoreGreen;
        bool mIgnoreBlue;
        
        // horizontal scaling
        double mLowerScale;
        double mUpperScale;
        
        // shortcut
        //
        // Given value X between mLowerBound and mUpperBound we need to find a value 
        // Z between mLowerScale and mUpperScale. This involves a constant multiplicaation
        // by (mUpperScale - mLowerScale)/(mUpperBound - mLowerBound). So we work this out once.
        //
        double mScaleMultiplier;
        double mTickSize;          // size of a tick
};

#endif //Rainbow_h

