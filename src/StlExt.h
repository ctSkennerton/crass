//
//  StlExt.h
//  crass
//
//  Created by Connor Skennerton on 1/07/11.
//  Copyright 2011 Australian Centre for Ecogenomics. All rights reserved.
//

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