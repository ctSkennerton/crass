//
//  DRArray.h
//  crass
//
//  Created by Connor Skennerton on 29/10/13.
//  Copyright (c) 2013 Australian Centre for Ecogenomics. All rights reserved.
//

#ifndef __crass__DRArray__
#define __crass__DRArray__
#include <vector>
#include <iostream>


// swap the integers without a tmp value
//http://graphics.stanford.edu/~seander/bithacks.html#SwappingValuesXOR
#define SWAP(a, b) (((a) ^ (b)) && ((b) ^= (a) ^= (b), (a) ^= (b)))

namespace crass {
        
    class RepeatArray {
        typedef std::vector<int> storage_t;

    public:
        class RepeatIterator {
        public:
            typedef RepeatIterator self_t;
            typedef int value_t;
            typedef int& reference_t;
            typedef int* pointer_t;
            typedef std::pair<reference_t, reference_t> return_t;
            
            RepeatIterator(){}
            
            RepeatIterator(storage_t::iterator input) : current_pos(input){}
            
            return_t operator *() {
                return return_t(*(current_pos), *(current_pos + 1 ));
            }
            
            self_t& operator++() { current_pos+=2; return *this; }
            self_t operator++(int ) { self_t i = *this; operator++(); return i; }
            
            self_t& operator--() { current_pos-=2; return *this; }
            self_t operator--(int) { self_t i = *this; operator--(); return i; }

            
            bool operator==(const self_t& rhs) { return current_pos == rhs.current_pos; }
            bool operator!=(const self_t& rhs) { return current_pos != rhs.current_pos; }
            bool operator<(const self_t& rhs) { return current_pos < rhs.current_pos; }
            bool operator<=(const self_t& rhs) { return current_pos <= rhs.current_pos; }
            bool operator>(const self_t& rhs) { return current_pos > rhs.current_pos; }
            bool operator>=(const self_t& rhs) { return current_pos >= rhs.current_pos; }
            
                
                self_t& operator+=(const value_t& offset) {
                    current_pos += (offset*2);
                    return *this;
                }
                self_t operator+(const value_t& offset) {
                    self_t ret = *this;
                    return ret += offset;
                }
                
                self_t& operator-=(const value_t& offset) {
                    current_pos -= (offset*2);
                    return *this;
                }
                self_t operator-(const value_t& offset) {
                    self_t ret = *this;
                    return ret -= offset;
                }

                
        private:
            storage_t::iterator current_pos;
        };
            
        class SpacerIterator {
        public:
            typedef SpacerIterator self_t;
            typedef int value_t;
            typedef int& reference_t;
            typedef int* pointer_t;
            typedef std::pair<reference_t, reference_t> return_t;
            
            SpacerIterator(){}
            
            SpacerIterator(storage_t::iterator input) : current_pos(input){}
            
            return_t operator *() {
                return return_t(*(current_pos), *(current_pos + 1 ));
            }
            
            self_t& operator++() { current_pos+=2; return *this; }
            self_t operator++(int ) { self_t i = *this; operator++(); return i; }

            self_t& operator--() { current_pos-=2; return *this; }
            self_t operator--(int) { self_t i = *this; operator--(); return i; }

            
            bool operator==(const self_t& rhs) { return current_pos == rhs.current_pos; }
            bool operator!=(const self_t& rhs) { return current_pos != rhs.current_pos; }
            bool operator<(const self_t& rhs) { return current_pos < rhs.current_pos; }
            bool operator<=(const self_t& rhs) { return current_pos <= rhs.current_pos; }
            bool operator>(const self_t& rhs) { return current_pos > rhs.current_pos; }
            bool operator>=(const self_t& rhs) { return current_pos >= rhs.current_pos; }
            
                
            self_t& operator+=(const value_t& offset) {
                current_pos += (offset*2);
                return *this;
            }
            self_t operator+(const value_t& offset) {
                self_t ret = *this;
                return ret += offset;
            }

            self_t& operator-=(const value_t& offset) {
                current_pos -= (offset*2);
                return *this;
            }
            self_t operator-(const value_t& offset) {
                self_t ret = *this;
                return ret -= offset;
            }

            
            
        private:
            storage_t::iterator current_pos;
        };
        
        /** add
         *  add in the position of a new repeat using 0-indexed half open addressing
         *
         *  @param start 0-indexed starting position of the repeat in the read
         *  @param end 0-indexed position pointing to 1 passed the end of the repeat  
         */
        void add(int start, int end) {
            positions.push_back(start);
            positions.push_back(end);
        }
        
        void dump() {
            for (auto i : positions) {
                std::cout <<i<<",";
            }
            std::cout <<std::endl;
        }
        
            bool empty() { return positions.empty();}
        
        
        int numberOfRepeats();
        int numberOfSpacers();
            
        void reverseRepeatPositions(int finalIndex);
        
        RepeatIterator repeatBegin(){return RepeatIterator(positions.begin());}
        RepeatIterator repeatEnd(){return RepeatIterator(positions.end());}
            
        SpacerIterator spacerBegin(){return SpacerIterator(positions.begin() + 1);}
        SpacerIterator spacerEnd(){return SpacerIterator(positions.end() - 1);}
            
        RepeatIterator repeatAt(int i);
        SpacerIterator spacerAt(int i);
            int repeatLengthAt(int i);
            int spacerLengthAt(int i);
            

    protected:
        std::vector<int> positions;
    };
    
    std::ostream& operator<<(std::ostream& output, RepeatArray& data);
}
#endif /* defined(__crass__DRArray__) */

