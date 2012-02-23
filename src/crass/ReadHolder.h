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

// typedefs
typedef std::vector<unsigned int> StartStopList;
typedef std::vector<unsigned int>::iterator StartStopListIterator;
typedef std::vector<unsigned int>::reverse_iterator StartStopListRIterator;




class ReadHolder 
{
    public:


        ReadHolder() 
        { 
            RH_LastDREnd = 0; 
            RH_NextSpacerStart = 0; 
            RH_isSqueezed = false;
            RH_IsFasta = true;
        }  
        
        ReadHolder(std::string s, std::string h) 
        {
            RH_Seq = s; 
            RH_Header = h; 
            RH_LastDREnd = 0; 
            RH_NextSpacerStart = 0; 
            RH_isSqueezed = false;
            RH_IsFasta = true;

        }

        ReadHolder(const char * s, const char * h) 
        {
            RH_Seq = s; 
            RH_Header = h;
            RH_LastDREnd = 0; 
            RH_NextSpacerStart = 0; 
            RH_isSqueezed = false;
            RH_IsFasta = true;

        }
        ReadHolder(std::string s, std::string h, std::string c, std::string q) 
        {
            RH_Seq = s; 
            RH_Header = h; 
            RH_Comment = c;
            RH_Qual = q;
            RH_LastDREnd = 0; 
            RH_NextSpacerStart = 0; 
            RH_isSqueezed = false;
            RH_IsFasta = false;
        }
        
        ReadHolder(const char * s, const char * h, const char * c, const char * q) 
        {
            RH_Seq = s; 
            RH_Header = h;
            RH_Comment = c;
            RH_Qual = q;
            RH_LastDREnd = 0; 
            RH_NextSpacerStart = 0; 
            RH_isSqueezed = false;
            RH_IsFasta = false;

        }
        
        ~ReadHolder(void)
        {
            clear();
        }
        
        void clear(void)
        {
            RH_Seq.clear();
            RH_StartStops.clear();
            RH_Header.clear();
            RH_Rle.clear();
            RH_Comment.clear();
            RH_Qual.clear();
        }


    
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
        inline bool getIsFasta(void)
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
        
        inline bool getSqueezed(void)
        {
            return this->RH_isSqueezed;
        }
        
        inline StartStopList getStartStopList(void)
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
    
        inline int start()
        {
            return this->front();
        }
        
        inline int getFirstRepeatStart()
        {
            return this->front();
        }
        
        // return the second last element in the start stop list
        // which is equal to the start of the last repeat
        inline int getLastRepeatStart()
        {
            StartStopListIterator iter = RH_StartStops.end() - 2;
            return *iter;
        }
        
        inline unsigned int numRepeats()
        {
            return (unsigned int)(RH_StartStops.size()/2);
        }
        
        inline unsigned int numSpacers()
        {
            return numRepeats() - 1;
        }
        
        unsigned int front(void)
        {
            return RH_StartStops.front();
        }
        unsigned int back(void)
        {
            return RH_StartStops.back();
        }
        
        int getSeqLength(void)
        {
            return (int)RH_Seq.length();
        }
        
        unsigned int getStartStopListSize(void)
        {
            return (int)RH_StartStops.size();
        }
    
        inline unsigned int getRepeatLength()
        {
            return RH_RepeatLength;
        }
    
        inline char getSeqCharAt(int i)
        {
            return RH_Seq[i];
        }
        
        unsigned int getRepeatAt(unsigned int i);

        std::string repeatStringAt(unsigned int i);
        
        std::string spacerStringAt(unsigned int i);
    
        unsigned int getAverageSpacerLength(void);
        
        std::vector<std::string> getAllSpacerStrings(void);
    
        std::vector<std::string> getAllRepeatStrings(void);
    
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
        inline void setSequence(std::string _sequence)
        {
            RH_RepeatLength = 0;
            RH_Seq = _sequence;
        }

        inline void setRepeatLength(int length)
        {
            RH_RepeatLength = length;
        }
        
        inline void incrementRepeatLength(void)
        {
            RH_RepeatLength++;
        }
    
        inline void decrementRepeatLength(void)
        {
            RH_RepeatLength--;
        }
    
        inline void setHeader(std::string h)
        {
            this->RH_Header = h;
        }
        
        inline void setDRLowLexi(bool b)
        {
            this->RH_WasLowLexi = b;
        }
        
        inline void setSqueezed(bool b)
        {
            this->RH_isSqueezed = b;
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


        int startStopsAt(int i)
        {
            return this->RH_StartStops.at(i);
        }        
        void startStopsAdd(unsigned int, unsigned int);


		void clearStartStops(void)
		{
			RH_StartStops.clear();
		}
		
		void removeRepeat(unsigned int val);
		
		void dropPartials(void);
	
		void reverseStartStops(void);           // fix start stops what got corrupted during revcomping
	
		// update the DR after finding the TRUE DR
		void updateStartStops(int frontOffset, std::string * DR, const options * opts);
		
		// the positions are the start positions of the direct repeats
		// 
		inline int repeatSpacing(unsigned int pos1, unsigned int pos2)
		{
			return (getRepeatAt(pos2) - getRepeatAt(pos1));
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
                
        std::string expand(void);               // returns a copy of the string with homopolymers
        
        std::string expand(bool fixStopStarts); // returns a copy of the string with homopolymers (fixes stop starts)
        
        // cut DRs and Specers
        bool getFirstDR(std::string * retStr);
        
        bool getNextDR(std::string * retStr);
        
        bool getFirstSpacer(std::string * retStr);
        
        bool getNextSpacer(std::string * retStr);        
    
        std::string toStringInColumns(void);    
    
        std::string splitApart(void);
        
        std::string splitApartSimple(void);
        
        void printContents(void);
        
        void logContents(int logLevel);
    
        inline std::ostream& print(std::ostream& s);
    
    private:
        // members
        std::string RH_Rle;                     // Run length encoded string
        std::string RH_Header;                  // Header for the sequence
        std::string RH_Comment;                 // The comment attribute of the sequence
        std::string RH_Qual;                    // The quality of the sequence
        bool RH_IsFasta;                        // boolean to tell us if the read is fastq or fasta
        std::string RH_Seq;                     // The DR_lowlexi sequence of this read
        bool RH_WasLowLexi;                     // was the sequence DR_low lexi in the file?
        StartStopList RH_StartStops;            // start stops for DRs, (must be even in length!)
        bool RH_isSqueezed;                     // Bool to tell whether the read has homopolymers removed
        int RH_LastDREnd;                       // the end of the last DR cut (offset of the iterator)
        int RH_NextSpacerStart;                 // the end of the last spacer cut (offset of the iterator)
        int RH_RepeatLength;
};

// overloaded operators 
std::ostream& operator<< (std::ostream& s,  ReadHolder& c);
#endif //ReadHolder_h
