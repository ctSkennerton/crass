//
//  SeqUtils.cpp
//  MCD
//
//  Created by Connor Skennerton on 31/05/11.
//  Copyright 2011 The Faculty of EAIT. All rights reserved.
//

#include "SeqUtils.h"
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;
string reverseComplement(string str)
{
    string revcomp_string = "";
    string::reverse_iterator rit;
    cout<<str<<endl;
    for ( rit = str.rbegin() ; rit < str.rend(); rit++ )
    {
		switch ((*rit)) 
        {
            case 'A':
            case 'a':
                revcomp_string += 'T';
                break;
            case 'C':
            case 'c':
                revcomp_string += 'G';
                break;
            case 'G':
            case 'g':
                revcomp_string += 'C';
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                revcomp_string += 'A';
                break;
            case 'M':
            case 'm':
                revcomp_string +='K'; 
                break;
            case 'R':
            case 'r':
                revcomp_string += 'Y';
                break;
            case 'W':
            case 'w':
                revcomp_string += 'W';
                break;
            case 'S':
            case 's':
                revcomp_string += 'S';
                break;
            case 'Y':
            case 'y':
                revcomp_string += 'R';
                break;
            case 'K':
            case 'k':
                revcomp_string += 'M';
                break;
            case 'V':
            case 'v':
                revcomp_string += 'B';
                break;
            case 'H':
            case 'h':
                revcomp_string += 'D';
                break;
            case 'D':
            case 'd':
                revcomp_string += 'H';
                break;
            case 'B':
            case 'b':
                revcomp_string += 'V';
                break;
            case 'N':
            case 'n':
                revcomp_string += 'N';
                break;
            default:
                revcomp_string += 'N';
                break;
		}
	}
    cout<<revcomp_string<<endl;
    return revcomp_string;
}

string laurenize (string seq1)
{
    string seq2 = reverseComplement(seq1);
    if (lexicographical_compare(seq1.begin(), seq1.end(), seq2.begin(), seq2.end()))
    {
        return seq1;
    }
    return seq2;
}