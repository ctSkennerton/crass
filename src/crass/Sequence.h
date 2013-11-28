//
//  Sequence.h
//  crass
//
//  Created by Connor Skennerton on 1/11/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#ifndef __crass__Sequence__
#define __crass__Sequence__

#include <iostream>
#include <string>

#include "RepeatArray.h"
namespace crass {

// global functions
std::string reverseComplement(std::string str);

std::string laurenize(std::string seq);
    
    class RawRead {

    public:
        
        class SpacerStringIterator {
        public:
            typedef SpacerStringIterator self_t;
            typedef std::string value_t;
            typedef value_t& reference_t;
            typedef value_t* pointer_t;
            
            SpacerStringIterator(std::string& input, RepeatArray::SpacerIterator si) : data(input), current_pos(si){}
            
            value_t operator *() {
                return data.substr((*current_pos).first, (*current_pos).second - (*current_pos).first);
            }
            
            self_t operator++() { self_t i = *this; ++current_pos; return i; }
            self_t operator++(int junk) { ++current_pos; return *this; }

            
            bool operator==(const self_t& rhs) { return current_pos == rhs.current_pos; }
            bool operator!=(const self_t& rhs) { return current_pos != rhs.current_pos; }
            bool operator<(const self_t& rhs) { return current_pos < rhs.current_pos; }
            bool operator<=(const self_t& rhs) { return current_pos <= rhs.current_pos; }
            bool operator>(const self_t& rhs) { return current_pos > rhs.current_pos; }
            bool operator>=(const self_t& rhs) { return current_pos >= rhs.current_pos; }
            
            
        private:
            std::string& data;
            RepeatArray::SpacerIterator current_pos;
        };
            
        class RepeatStringIterator {
        public:
            typedef RepeatStringIterator self_t;
            typedef std::string value_t;
            typedef value_t& reference_t;
            typedef value_t* pointer_t;
            
            RepeatStringIterator(std::string& input, RepeatArray::RepeatIterator si) : data(input), current_pos(si){}
            
            value_t operator *() {
                return data.substr((*current_pos).first, (*current_pos).second - (*current_pos).first);
            }
            
            self_t operator++() { self_t i = *this; ++current_pos; return i; }
            self_t operator++(int junk) { ++current_pos; return *this; }
            
            
            bool operator==(const self_t& rhs) { return current_pos == rhs.current_pos; }
            bool operator!=(const self_t& rhs) { return current_pos != rhs.current_pos; }
            bool operator<(const self_t& rhs) { return current_pos < rhs.current_pos; }
            bool operator<=(const self_t& rhs) { return current_pos <= rhs.current_pos; }
            bool operator>(const self_t& rhs) { return current_pos > rhs.current_pos; }
            bool operator>=(const self_t& rhs) { return current_pos >= rhs.current_pos; }
                
                
            private:
                std::string& data;
                RepeatArray::RepeatIterator current_pos;
            };
        
        RawRead()
        {}
        
        RawRead(std::string n, std::string c, std::string s, std::string q, RepeatArray r) :
                mName(n),
                mComment(c),
                mSeq(s),
                mQual(q),
                mRepeatPositions(r),
                mRepeatLowLexi(false),
                mPrevRead(nullptr),
                mNextRead(nullptr)
        {}
        
        ~RawRead()
        {}
        
        // get/set
        void identifier(std::string i) {
            mName = i;
        }
        std::string identifier() {
            return mName;
        }
        
        void comment(std::string i) {
            mComment = i;
        }
        std::string comment() {
            return mComment;
        }
        
        void seq(std::string i) {
            mSeq = i;
        }
        std::string seq() {
            return mSeq;
        }
        
        void quality(std::string i) {
            mQual = i;
        }
        std::string quality() {
            return mQual;
        }
        
        void repeatLowestLexicographicalForm(bool b) {
            mRepeatLowLexi = b;
        }
        bool repeatLowestLexicographicalForm() {
            return mRepeatLowLexi;
        }
            
        void nextRead(RawRead *n) {
            mNextRead = n;
        }
        RawRead * nextRead() {
            return mNextRead;
        }
        
        void previousRead(RawRead *p) {
            mPrevRead = p;
        }
        RawRead * previousRead() {
            return mPrevRead;
        }
        
        void push_back(int s, int e) {
            mRepeatPositions.add(s, e);
        }
            
        bool firstReadInFragment() {
            return (mNextRead != nullptr) && (mPrevRead == nullptr);
        }
        
        bool lastReadInFragment() {
            return (mPrevRead != nullptr) && (mNextRead == nullptr);
        }
        
        void inspect() {
            std::cout<<"mName="<<mName<<"\n";
            std::cout<<"mComment="<<mComment<<"\n";
            std::cout<<"mSeq="<<mSeq<<"\n";
            std::cout<<"mQual="<<mQual<<"\n";
            std::cout<<"mRepeatLowLexi="<<mRepeatLowLexi<<"\n";
            std::cout<<mRepeatPositions<<std::endl;
            
        }
    
        size_t length() {
            return mSeq.length();
        }
        size_t size() {
            return mSeq.length();
        }
            
            void clearRepeatPositions() {mRepeatPositions.clear();}
            
        int numberOfRepeats() {
            return mRepeatPositions.numberOfRepeats();
        }
        
        int numberOfSpacers() {
            return mRepeatPositions.numberOfSpacers();
        }
            
        bool startPartial() {
            auto it = mRepeatPositions.repeatBegin();
            return (*it).first == 0;
        }
        
        bool endPartial() {
            auto it = mRepeatPositions.repeatEnd() - 1;
            return (*it).first == length();
        }
            
        
        RepeatArray::RepeatIterator repeatBegin(){return mRepeatPositions.repeatBegin();}
        RepeatArray::RepeatIterator repeatEnd(){return mRepeatPositions.repeatEnd();}
        
        RepeatArray::SpacerIterator spacerBegin(){return mRepeatPositions.spacerBegin();}
        RepeatArray::SpacerIterator spacerEnd(){return mRepeatPositions.spacerEnd();}
            
            RepeatArray::RepeatIterator repeatAt(int i){return mRepeatPositions.repeatAt(i);}
            RepeatArray::SpacerIterator spacerAt(int i){return mRepeatPositions.spacerAt(i);}
        
        RawRead::SpacerStringIterator spacerStringBegin(){return SpacerStringIterator(mSeq, mRepeatPositions.spacerBegin());}
        RawRead::SpacerStringIterator spacerStringEnd(){return SpacerStringIterator(mSeq, mRepeatPositions.spacerEnd());}
        
        RawRead::RepeatStringIterator repeatStringBegin(){return RepeatStringIterator(mSeq, mRepeatPositions.repeatBegin());}
        RawRead::RepeatStringIterator repeatStringEnd(){return RepeatStringIterator(mSeq, mRepeatPositions.repeatEnd());}
        
        RawRead::RepeatStringIterator repeatStringAt(int i){return RepeatStringIterator(mSeq, mRepeatPositions.repeatAt(i));}
        RawRead::SpacerStringIterator spacerStringAt(int i){return SpacerStringIterator(mSeq, mRepeatPositions.spacerAt(i));}
        
        char& operator[](int i) {
            return mSeq[i];
        }
    
        void revComp();
            
        std::string orientateRepeatLowLexi();
            
        int getFirstNonPartialRepeatLength();
            
        RawRead::RepeatStringIterator getFirstNonPartialRepeat();
            RepeatArray::RepeatIterator getFirstNonPartialRepeatPositions();
            



        
    private:
        
        //members
        std::string mName;
        std::string mComment;
        std::string mSeq;
        std::string mQual;
        RepeatArray mRepeatPositions;
        bool mRepeatLowLexi;
        RawRead *mPrevRead;
        RawRead *mNextRead;
    };

}
#endif /* defined(__crass__RawRead__) */
