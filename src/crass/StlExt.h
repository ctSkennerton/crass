/*
 * StlExt.h
 *
 * Copyright (C) 2011 - Connor Skennerton
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef STLEXT_H
#define STLEXT_H

#include <sstream>
#include <map>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <vector>
#include <iterator>
#include <numeric>
#include <iostream>



template <class T1, class T2>
void addOrIncrement(std::map<T1, T2> &inMap, T1 &searchThing)
{
    
    if (inMap.find(searchThing) != inMap.end())
    {
        inMap[searchThing] += 1;
    }
    else
    {
        inMap[searchThing] = 1;
    }
}

template <class T>
inline std::string to_string (const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

template <class T1, class T2, class T3>
void mapToVector(std::map<T1, T2>& map, std::vector<T3>& vector) 
{
    
    typename std::map<T1, T2>::iterator iter = map.begin();
    while (iter != map.end()) 
    {
        vector.push_back(iter->first);
        iter++;
    }
}

// templated function to split a string on delimeters and return it in a container of class T
template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", const bool trimEmpty = false)
{
    std::string::size_type pos, lastPos = 0;
    while(true)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos)
        {
            pos = str.length();
            
            if(pos != lastPos || !trimEmpty)
                tokens.push_back( typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos - lastPos ));
            
            break;
        }
        else
        {
            if(pos != lastPos || !trimEmpty)
                tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos ));
        }
        
        lastPos = pos + 1;
    }
}

template < class ContainerT >
void split(const std::string& str, ContainerT& tokens, const std::string& delimiters = ",", const bool trimEmpty = false)
{
    std::string::size_type pos, lastPos = 0;
    while(true)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos)
        {
            pos = str.length();
            
            if(pos != lastPos || !trimEmpty)
                tokens.insert( typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos - lastPos ));
            
            break;
        }
        else
        {
            if(pos != lastPos || !trimEmpty)
                tokens.insert(typename ContainerT::value_type(str.data()+lastPos, (typename ContainerT::value_type::size_type)pos-lastPos ));
        }
        
        lastPos = pos + 1;
    }
}

// class T1 is the container that holds the numbers and class T2 is the return type like int or double
template <class T >
typename T::value_type mean(T& container) {
    return (std::accumulate(container.begin(), container.end(), static_cast<typename T::value_type>(0))/container.size());
}

template <class T>
typename T::value_type median(T& container) {
    std::nth_element(container.begin(), container.begin()+container.size()/2, container.end());
    return *(container.begin()+container.size()/2);
}

template <class T>
typename T::value_type percentile(T& container, double percentile) {
    std::nth_element(container.begin(), container.begin()+container.size()*percentile, container.end());
    return *(container.begin()+container.size()*percentile);
}

template <class T>
typename T::value_type mode(T& container) {
    std::vector<typename T::value_type> histogram(container.size(),0);
    typename T::iterator iter = container.begin();
    while (iter != container.end()) {
        histogram[*iter++]++;
    }
    return std::max_element(histogram.begin(), histogram.end()) - histogram.begin();
}

template <class T>
double standardDeviation(T& container) {
    double average = static_cast<double>( mean(container));
    std::vector<double> temp;
    typename T::iterator iter;
    for (iter = container.begin(); iter != container.end(); iter++) {
        double i = static_cast<double>(*iter) - average;
        temp.push_back(i*i);
    }
    return std::sqrt(std::accumulate(temp.begin(), temp.end(), static_cast<double>(0))/temp.size() );
}

#endif
