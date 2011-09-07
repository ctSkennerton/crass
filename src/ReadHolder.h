// File: ReadHolder.h
// Original Author: Michael Imelfort 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
// Holder of reads. Identified by the various serach algorithms
// Storage class, so stupidly public!
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

#ifndef ReadHolder_h
#define ReadHolder_h

// system includes
#include <iostream>
#include <vector>
#include <map>
// local includes
#include "crassDefines.h"

// typedefs
typedef std::vector<unsigned int> StartStopList;
typedef std::vector<unsigned int>::iterator StartStopListIterator;
typedef std::vector<unsigned int>::reverse_iterator StartStopListRIterator;


class ReadHolder 
{
    public:
        
        ReadHolder() 
        { 
            mLastDREnd = 0; 
            mLastSpacerEnd = 0; 
            RH_isSqueezed = false;
        }  
        
        ReadHolder(std::string& s, std::string& h) 
        {
            RH_Seq = s; 
            RH_Header = h; 
            mLastDREnd = 0; 
            mLastSpacerEnd = 0; 
            RH_isSqueezed = false;
        }

        ReadHolder(const char * s, const char * h) 
        {
            RH_Seq = s; 
            RH_Header = h;
            mLastDREnd = 0; 
            mLastSpacerEnd = 0; 
            RH_isSqueezed = false;
        }
        
        void clear(void)
        {
            RH_Seq.clear();
            RH_StartStops.clear();
            RH_Header.clear();
        }
        
        //----
        // Getters
        //
        std::string seq(void)
        {
            return this->RH_Seq;
        }
    
        std::string header(void)
        {
            return this->RH_Header;
        }
        
        std::string seqLiteral(void)
        {
            return this->RH_Seq;
        }
    
        bool isLowLexi(void)
        {
            return this->RH_WasLowLexi;
        }
        
        bool isSqueezed(void)
        {
            return this->RH_isSqueezed;
        }
        
        StartStopList drPos(void)
        {
            return this->RH_StartStops;
        }
        
        int lastDRPos(void)
        {
            return this->mLastDREnd;
        }
        
        int lastSpacerPos(void)
        {
            return this->mLastSpacerEnd;
        }
        
        //----
        //setters
        // 
        void seq(std::string s)
        {
            this->RH_Seq = s;
        }

        void header(std::string h)
        {
            this->RH_Header = h;
        }
        
        void isLowLexi(bool b)
        {
            this->RH_WasLowLexi = b;
        }
        
        void isSqueezed(bool b)
        {
            this->RH_isSqueezed = b;
        }
        
        void drPos(StartStopList& l)
        {
            this->RH_StartStops.clear();
            this->RH_StartStops = l;
        }
        
        void lastDRPos(int i)
        {
            this->mLastDREnd = i;
        }
        
        void lastSpacerPos(int i)
        {
            this->mLastSpacerEnd = i;
        }
        
        int seqLength(void)
        {
            return (int)this->RH_Seq.length();
        }
        
        int drListSize(void)
        {
            return (int)this->RH_StartStops.size();
        }
        
        //----
        // Element access to the start stop list
        //
        int at(int i)
        {
            return this->RH_StartStops.at(i);
        }        
        void add(int, int);
        
        int front(void)
        {
            return this->RH_StartStops.front();
        }
        int back(void)
        {
            return this->RH_StartStops.back();
        }
    
        StartStopListIterator begin(void)
        {
            return this->RH_StartStops.begin();
        }
        
        StartStopListIterator end(void)
        {
            return this->RH_StartStops.end();
        }
    
        unsigned int& operator[]( const unsigned int i)
        {
            return this->RH_StartStops[i];
        }
    
        void reverseComplementSeq(void);        // reverse complement the sequence and fix the start stops
        void reverseStartStops(void);           // fix start stops what got corrupted during revcomping
        
        
        void encode(void);                      // transforms a sequence to a run length encoding form
        void decode (void);                     // transforms a sequence back to its original state
        std::string squeeze(void);              //returns a copy of the string without homopolymers
        std::string expand(void);               // returns a copy of the string with homopolymers
        std::string expand(bool fixStopStarts); // returns a copy of the string with homopolymers (fixes stop starts)
        
        // update the DR after finding the TRUE DR
        void updateStartStops(int frontOffset, std::string * DR, const options * opts);

        // cut DRs and Specers
        bool getFirstDR(std::string * retStr);
        bool getNextDR(std::string * retStr);
        bool getFirstSpacer(std::string * retStr);
        bool getNextSpacer(std::string * retStr);        
        std::string splitApart(void);
        std::string splitApartSimple(void);
        
        void printContents(void);
        void logContents(int logLevel);
    private:
        // members
        std::string RH_Rle;                     // Run length encoded string
        std::string RH_Header;                  // Header for the sequence
        std::string RH_Seq;                     // The DR_lowlexi sequence of this read
        bool RH_WasLowLexi;                     // was the sequence DR_low lexi in the file?
        StartStopList RH_StartStops;            // start stops for DRs, (must be even in length!)
        bool RH_isSqueezed;                     // Bool to tell whether the read has homopolymers removed
        int mLastDREnd;                         // the end of the last DR cut (offset of the iterator)
        int mLastSpacerEnd;                     // the end of the last spacer cut (offset of the iterator)
};

#endif //ReadHolder_h
