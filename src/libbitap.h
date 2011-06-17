/*
 *  bitap.h is part of the MCD project
 *  
 *  Created by Connor Skennerton.
 *  Copyright 2011 Connor Skennerton & Michael Imelfort. All rights reserved. 
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

#ifndef LIBBITAP_H
#define LIBBITAP_H
#define LIBBITAP_MAJOR_VERSION 1
#define LIBBITAP_MINOR_VERSION 2

//#include <string.h>

typedef enum {
    anyJump, infixLeft, infixRight
} jumpType;

struct bitapJumpType {
    int foffset, fmask, toffset, tmask;
    jumpType t;
};
/* This structure contains data that is private to libbitap and should not */
/* be modified directly */
typedef struct {
    int l, n;            /* l is the number of states. n is the number of */
    /* 'or'ing operations stored. */
    int *s, *x;          /* s shows which letter is allowed in which place. */
    /* x determines in which areas exact matches are */
    struct bitapJumpType *j; /* required. */
    int approxMode;          /* Only used during parsing. To determine x */
} bitapType;

int NewBitap (bitapType *b, const char *regex);

const char * FindWithBitap (bitapType *b, const char *in, int inl, int e, int *ereturn, const char **bReturn);

void DeleteBitap (bitapType *b);

#endif
