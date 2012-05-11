/*
 *  StatsManager.h is part of the crass project
 *  
 *  Created by Connor Skennerton on 8/01/12.
 *  Copyright 2012 Connor Skennerton. All rights reserved. 
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

#ifndef crass_StatsManager_h
#define crass_StatsManager_h
/*
 This is a simple little class for calculating the mean etc... 
 for any container that accepts a push_back funtion (vector, deque, list, stack, queue)
 */
#include <cmath>
#include <numeric>
#include <algorithm>

template <class T>
class StatsManager {
    T container;

public:
    StatsManager (void){}
    typename T::value_type mean(void);
    typename T::value_type median(void);
    typename T::value_type percentile(double percent);
    typename T::value_type mode(void);
    double standardDeviation(void);
    void add(typename T::value_type a);
    bool remove(typename T::value_type a);
    void clear(void);
    
};

template <class T >
typename T::value_type StatsManager<T>::mean() {
    return (std::accumulate(container.begin(), container.end(), static_cast<typename T::value_type>(0))/container.size());
}

template <class T>
typename T::value_type StatsManager<T>::median(void) {
    std::nth_element(container.begin(), container.begin()+container.size()/2, container.end());
    return *(container.begin()+container.size()/2);
}

template <class T>
typename T::value_type StatsManager<T>::percentile( double percentile) {
    std::nth_element(container.begin(), container.begin()+container.size()*percentile, container.end());
    return *(container.begin()+container.size()*percentile);
}

template <class T>
typename T::value_type StatsManager<T>::mode(void) {
    std::vector<typename T::value_type> histogram(container.size(),0);
    typename T::iterator iter = container.begin();
    while (iter != container.end()) {
        histogram[*iter++]++;
    }
    return std::max_element(histogram.begin(), histogram.end()) - histogram.begin();
}

template <class T>
double StatsManager<T>::standardDeviation(void) {
    double average = static_cast<double>( mean());
    std::vector<double> temp;
    typename T::iterator iter;
    for (iter = container.begin(); iter != container.end(); iter++) {
        double i = static_cast<double>(*iter) - average;
        temp.push_back(i*i);
    }
    return std::sqrt(std::accumulate(temp.begin(), temp.end(), static_cast<double>(0))/temp.size() );

    
}

template <class T>
bool StatsManager<T>::remove(typename T::value_type a) {
    typename T::iterator iter;
    bool success = false;
    for (iter = container.begin(); iter != container.end(); iter++) {
        if (*iter == a) {
            container.erase(*iter);
            success = true;
            break;
        }
    }
    return success;
}

template <class T>
void StatsManager<T>::add(typename T::value_type a) {
    container.push_back(a);
}

template <class T>
void StatsManager<T>::clear(void) {
    container.clear();
}
#endif
