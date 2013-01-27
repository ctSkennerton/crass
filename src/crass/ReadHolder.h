// File: ReadHolder.h
// Original Author: Michael Imelfort 2011
// Hacked and Extended: Connor Skennerton 2011
// --------------------------------------------------------------------
//
// OVERVIEW:
//
//  Main class for all things read related.  Holds information about the
//  sequence and the direct repeats contained in the read
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
#include "DRArray.h"

// typedefs
typedef crispr::RepeatArray<unsigned int> StartStopList;
typedef StartStopList::iterator StartStopListIterator;
//typedef StartStopList::reverse_iterator StartStopListRIterator;

namespace crass {
    class ReadHolder {
    public:
        
        ReadHolder(std::string s, std::string h): 
        RH_LastDREnd(0), 
        RH_NextSpacerStart(0), 
        RH_isSqueezed(false), 
        RH_IsFasta(true), 
        RH_Seq(s),
        RH_Header(h),
        RH_StartStops(s)
        {}
        
        ReadHolder(std::string s, std::string h, std::string c, std::string q): 
        RH_LastDREnd(0), 
        RH_NextSpacerStart(0), 
        RH_isSqueezed(false), 
        RH_IsFasta(false), 
        RH_Seq(s),
        RH_Header(h),
        RH_Comment(c),
        RH_Qual(q),
        RH_StartStops(s)
        {}
        
        //----
        // Getters
        //
        inline std::string getComment(void)
        {
            return this->RH_Comment;
        }
        inline std::string getQual(void)
        {
            return this->RH_Qual;
        }
        inline bool isFasta(void)
        {
            return this->RH_IsFasta;
        }
        inline std::string getSeq(void)
        {
            return this->RH_Seq;
        }
        
        inline std::string getHeader(void)
        {
            return this->RH_Header;
        }
        
        inline std::string getSeqRle(void)
        {
            return this->RH_Rle;
        }
        
        inline bool getLowLexi(void)
        {
            return this->RH_WasLowLexi;
        }
        
        inline bool isSqueezed(void)
        {
            return this->RH_isSqueezed;
        }
        
        inline crispr::RepeatArray<unsigned int> getStartStopList(void)
        {
            return this->RH_StartStops;
        }
        
        inline int getLastDRPos(void)
        {
            return this->RH_LastDREnd;
        }
        
        inline int getLastSpacerPos(void)
        {
            return this->RH_NextSpacerStart;
        }
        
        inline unsigned int start()
        {
            return this->front().first;
        }
        
        inline unsigned int getFirstRepeatStart()
        {
            return this->front().first;
        }
        
        // return the second last element in the start stop list
        // which is equal to the start of the last repeat
        inline unsigned int getLastRepeatStart()
        {
            return this->back().first;
        }
        
        inline unsigned int numRepeats()
        {
            return static_cast<unsigned int>(RH_StartStops.size());
        }
        
        inline unsigned int numSpacers()
        {
            return static_cast<unsigned int>(RH_StartStops.numberOfSpacers());
        }
        
        std::pair<unsigned int, unsigned int> front(void)
        {
            return RH_StartStops.front();
        }
        
        std::pair<unsigned int, unsigned int> back(void)
        {
            return RH_StartStops.back();
        }
        
        int getSeqLength(void)
        {
            return (int)RH_Seq.length();
        }
        
        unsigned int getStartStopListSize(void)
        {
            return static_cast<unsigned int>(RH_StartStops.size());
        }
        
        inline unsigned int getRepeatLength()
        {
            return RH_RepeatLength;
        }
        
        inline char getSeqCharAt(int i)
        {
            return RH_Seq[i];
        }
        
        inline std::pair<unsigned int, unsigned int>& getRepeatAt(unsigned int i) {
            return RH_StartStops.at(i);
        }
        
        inline std::string repeatStringAt(unsigned int i) {
            return RH_StartStops.repeatStringAt(i);
        }
        
        inline std::string spacerStringAt(unsigned int i) {
            return RH_StartStops.spacerStringAt(i);
        }
        
        unsigned int getAverageSpacerLength(void);
        
        void getAllSpacerStrings(std::vector<std::string>& spacers);
        
        void getAllRepeatStrings(std::vector<std::string>& repeats);
        
        int averageRepeatLength(void);

        //----
        //setters
        // 
        inline void setComment(std::string _comment)
        {
            RH_Comment = _comment;
        }
        inline void setQual(std::string _qual)
        {
            RH_Qual = _qual;
            RH_IsFasta = false;
        }

        inline void setRepeatLength(int length)
        {
            RH_RepeatLength = length;
        }
           
        inline void setHeader(std::string h)
        {
            this->RH_Header = h;
        }
        
        inline void setDRLowLexi(bool b)
        {
            this->RH_WasLowLexi = b;
        }
        
        void setLastDRPos(int i)
        {
            this->RH_LastDREnd = i;
        }
        
        void setLastSpacerPos(int i)
        {
            this->RH_NextSpacerStart = i;
        }
        //----
        // Element access to the start stop list
        //


        std::pair<unsigned int, unsigned int> startStopsAt(int i)
        {
            return this->RH_StartStops.at(i);
        }        
        void startStopsAdd(unsigned int, unsigned int);


        void clearStartStops(void)
        {
            RH_StartStops.clear();
        }
        
        inline void incrementRepeatLength(void)
        {
            RH_RepeatLength++;
        }
    
        inline void decrementRepeatLength(void)
        {
            RH_RepeatLength--;
        }

        void dropPartials(void);
    
        void reverseStartStops(void);           // fix start stops what got corrupted during revcomping
    
        // update the DR after finding the TRUE DR
        void updateStartStops(int frontOffset, std::string * DR, const options * opts);
        
        // the positions are the start positions of the direct repeats
        // 
        inline int repeatSpacing(unsigned int pos1, unsigned int pos2)
        {
            return (getRepeatAt(pos2).first - getRepeatAt(pos1).second);
        }

        crispr::RepeatArray<unsigned int>::iterator begin(void)
        {
            return this->RH_StartStops.begin();
        }
        
        crispr::RepeatArray<unsigned int>::iterator end(void)
        {
            return this->RH_StartStops.end();
        }


        std::pair<unsigned int, unsigned int>& operator[]( const unsigned int i)
        {
            return this->RH_StartStops[i];
        }


        std::string substr(int i, int j)
        {
            return RH_Seq.substr(i, j);
        }
        
        std::string substr(int i)
        {
            return RH_Seq.substr(i);
        }
        
        std::string substr(unsigned int i, unsigned int j)
        {
            return RH_Seq.substr(i, j );
        }
        
        std::string substr(unsigned int i)
        {
            return RH_Seq.substr(i);
        }
        
        std::string substr(size_t i, size_t j)
        {
            return RH_Seq.substr(i, j );
        }
        
        std::string substr(size_t i)
        {
            return RH_Seq.substr(i);
        }
    
        std::string DRLowLexi(void);            // Put the sequence in the form that makes the DR in it's laurenized form
    
        void reverseComplementSeq(void);        // reverse complement the sequence and fix the start stops
        
        void encode(void);                      // transforms a sequence to a run length encoding form
        
        void decode (void);                     // transforms a sequence back to its original state
                        
        std::string expand(bool fixStopStarts=false); // returns a copy of the string with homopolymers (fixes stop starts)
        
        // cut DRs and Specers
        bool getFirstDR(std::string * retStr);
        
        bool getNextDR(std::string * retStr);
        
        bool getFirstSpacer(std::string * retStr);
        
        bool getNextSpacer(std::string * retStr);        
    
        std::string toStringInColumns(void);    
    
        std::string splitApart(void);
        
        std::string splitApartSimple(void);
        
        void printContents(void);
        void printContents(std::ostream& out);
        
        void logContents(int logLevel);
    
        inline std::ostream& print(std::ostream& s);
    private:

        inline void setSqueezed(bool b)
        {
            this->RH_isSqueezed = b;
        }

        // members
        int RH_LastDREnd; 
        int RH_NextSpacerStart;
        bool RH_isSqueezed; 
        bool RH_IsFasta; 
        std::string RH_Seq;
        std::string RH_Header;
        std::string RH_Comment;                 // The comment attribute of the sequence
        std::string RH_Qual;                    // The quality of the sequence
        crispr::RepeatArray<unsigned int> RH_StartStops;
        std::string RH_Rle;                     // Run length encoded string
        bool RH_WasLowLexi;                     // was the sequence DR_low lexi in the file?
        int RH_RepeatLength;
        
    };
}

// overloaded operators 
std::ostream& operator<< (std::ostream& s,  crass::ReadHolder& c);

#endif //ReadHolder_h
