/*
 *  DRArray is part of the CRisprASSembler project
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

#ifndef crass_DRArray_h
#define crass_DRArray_h

#include <utility>
#include <vector>
#include <string>

namespace crispr {

    /** Storage class for individual repeats or spacers.
        This class can be used directly but it's most common
        usage will be as the return value from functions of 
        the RepeatArray class
     */
    template <class T>
    class ArraySubsequence {

    public:
        /** The type of CRISPR subunit.
         This enum should always be set to one of repeat or spacer.
         */
        enum unit_type{ REPEAT, SPACER };
        /**
         */
        ArraySubsequence(T s, T e, unit_type u, std::string& sequence) : 
        start(s), 
        stop(e), 
        type(u)
        {
            subseq = sequence.substr(start, stop - start + 1);
        }
        
        // member variables
        T start; /**< The start index of the sequence. This position is zero indexed and inclusive */
        T stop;  /**< The end index of the sequence. This position is zero indexed and inclusive */
        std::string subseq; /**< The DNA sequence corresponding to the start and stop positions. */  
        unit_type type; /**< The type of subsequence.  This should be one of repeat or spacer using the unit_type enum*/
    };
    
    
    template <class T>
    class RepeatArray {
        typedef std::pair<T, T> repeat;
        typedef std::vector< repeat > storage_type;
        typedef typename storage_type::iterator iterator_type;
        typedef typename storage_type::reverse_iterator reverse_iterator_type;
        
        // members
        storage_type container; 
        std::string& sequence;  // defined in a higher class like ReadHolder
    public:
        typedef typename storage_type::iterator iterator;

        
        RepeatArray(std::string& s) : sequence(s){}
        ~RepeatArray(){}
        
        
        inline typename std::vector<T>::size_type size() {
            return container.size();
        }
        
        inline typename std::vector<T>::size_type length() {
            return this->size();
        }
        
        inline typename std::vector<T>::size_type numberOfRepeats() {
            return this->size();
        }
        
        inline typename std::vector<T>::size_type numberOfSpacers() {
            return this->size() - 1;
        }
        
        //
        // element access
        //
        inline repeat front() {
            return container.front();
        }
        
        inline repeat back() {
            return container.back();
        }
        
        inline repeat& operator[](size_t pos) {
            return container[pos];
        }
        
        inline repeat& at(size_t pos) {
            return container.at(pos);
        }
        
        ArraySubsequence<T> repeatAt(size_t pos) {
            repeat tmp = this->at(pos);
            //std::string subseq(sequence, tmp.first, tmp.second - tmp.first + 1);
            ArraySubsequence<T> dr(tmp.first, tmp.second, ArraySubsequence<T>::REPEAT, sequence);
            return dr;
        }
        
        std::string repeatStringAt(size_t pos) {
            return this->repeatAt(pos).subseq;
        }
        
        ArraySubsequence<T> spacerAt(size_t pos) {
            repeat tmp = this->at(pos);
            repeat tmp2 = this->at(pos + 1);
            //std::string subseq(sequence, tmp.second + 1, tmp2.first - tmp.second + 1);
            ArraySubsequence<T> sp(tmp.second + 1, tmp2.first - 1, ArraySubsequence<T>::SPACER, sequence);
            return sp;
        }
        
        std::string spacerStringAt(size_t pos) {
            return this->spacerAt(pos).subseq;
        }
        
        //
        //  Element addition
        //
        inline void push_back(T& value1, T& value2 ) {
            container.push_back(std::make_pair(value1, value2));
        }
        
        inline void push_back(repeat value ) {
            container.push_back(value);
        }
        
        inline void insert(iterator_type pos, repeat value) {
            container.insert(pos, value);
        }
        //
        // element removal
        //
        void clear() {
            container.clear();
        }
        
        void erase(size_t pos) {
            iterator_type iter = this->begin() + pos;
            container.erase(iter);
        }
        
        void erase (iterator_type i) {
            container.erase(i);
        }
        
        void erase (iterator_type i, iterator_type j) {
            container.erase(i, j);
        }
        
        void reverse() {
            
            storage_type tmp;
            //StartStopList tmp_ss;
            T final_index = static_cast<T>(sequence.length() - 1);
            //int true_start_offset = (int)sequence.length() - this->back().second - 1;
            
            typename storage_type::reverse_iterator riter;
            for(riter = container.rbegin(); riter != container.rend(); riter++) {
                tmp.push_back(std::make_pair(final_index - riter->second, final_index - riter->first));
            }
            container.swap(tmp);            
        }
        
        //
        // Iterators
        //
        
        inline iterator_type begin() {
            return container.begin();
        }
        
        inline iterator_type end() {
            return container.end();
        }
        
        inline reverse_iterator_type rbegin() {
            return container.rbegin();
        }
        
        inline reverse_iterator_type rend() {
            return container.rend();
        }
        
    };
}


#endif
