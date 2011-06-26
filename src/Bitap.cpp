/*
 *  bitap.cpp is part of the CRisprASSembler project
 *  
 *  Created by Nic Roets http://rational.co.za/libbitap/
 *  Copyright 2005 Nic Roets. All rights reserved. 
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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#include "libbitap.h"

#define BITS (sizeof (int) * 8)
#define ALPHABET (1 << (sizeof (char) * 8))

static void NewJump (bitapType *b, int from, int to, jumpType t)
{
    if (b->s) {
        b->j[b->n].foffset = from / BITS;
        b->j[b->n].fmask = 1 << (from % BITS);
        b->j[b->n].toffset = to / BITS;
        b->j[b->n].tmask = 1 << (to % BITS);
        b->j[b->n].t = t;
    }
    b->n++;
}

static int ParseRegex (bitapType *b, const char *regex, const char *end)
{
    int lastLen = 0, lastL = b->l - 1, infixStart = b->l - 1, inverse, i;
    int infixEnd = -1; 
    /* Negative while no '|' was found, otherwise it marks the last bit before
     * the '|' which needs to be copied over to the end of the expression. */
    const char *start = regex;
    while ((!end || regex < end) && *regex != '\0' && *regex != ')') {
        if (*regex == '(') {
            lastL = b->l - 1;
            lastLen = ParseRegex (b, ++regex, end);
            if (lastLen < 0) return lastLen - (regex - start); /* Pass through */
            regex += lastLen;
            if (*regex++ != ')') return start - regex - 1;
            /* if the '(' was only matched by a '\0' */
            lastLen += 2; /* Account for the '(' and the ')' */
            continue;
        }
        if (*regex == '{') {
            for (i = 1; i < atoi (regex + 1); i++) {
                lastL = b->l;
                lastLen = ParseRegex (b, regex - lastLen, regex);
                /* This regex has already been scanned and will not contain errors */
            }
            
            for (regex++; regex < end && isdigit (*regex);) regex++;
            if (regex[strspn (regex, "0123456789")] == ',') {
                if (isdigit (regex[strspn (regex, "0123456789") + 1])) {
                    for (; i <= atoi (regex + strspn (regex, "0123456789") + 1); i++) {
                        lastL = b->l;
                        ParseRegex (b, regex - lastLen, regex); /* Already scanned */
                        NewJump (b, lastL - 1, b->l - 1, anyJump);
                    }
                }
                else NewJump (b, b->l - 1, lastL - 1, anyJump); /* The {xx,} case */
            }
            while (isdigit (*++regex) || *regex == ',') {}
            if (*regex++ != '}') return start - regex - 1;
            lastLen = 0; /* TODO: Handle 'a{1}*' */
            lastL = b->l - 1;
            continue;
        }
        if (*regex == '|') {
            regex++;
            if (infixEnd >= 0) NewJump (b, infixEnd - 1, b->l - 1, infixRight);
            /* If there are two or more '|' operators. End of left hand to end */
            
            infixEnd = b->l - 1;
            lastL = b->l; /* TODO: Handle 'abc|+' ? */
            NewJump (b, infixStart, b->l, infixLeft); /* beginning LH to RH */
            /* if (b->s) b->x[b->l / BITS] |= 1 << (b->l % BITS); */
            b->l++; /* Create an unmatchable char to stop shift from LHS to RHS */
            continue;
        }
        if (*regex == '+' || *regex == '*' || *regex == '?') 
        {
            if (b->l > 0) 
            { /* If regex start with '+' or '*', just ignore it */
                if (*regex != '?') NewJump (b, b->l - 1, lastL, anyJump);
                /* 'do{}' in C like '+' requires only 1 jump */
                if (*regex != '+') NewJump (b, lastL, b->l - 1, anyJump);
                /* 'if(){}' in C like '?' requires only 1 jump */
                /* 'while(){}' in C like '*' requires both jumps */
            }
            regex++;
            lastLen = 0; /* TODO : Handle a** */
            continue;
        }
        lastLen = 1;
        if (*regex == '[') {
            if ((inverse = *++regex == '^')) regex++;
            do {
                if (*regex == '\0') return start - regex - 1;
                i = regex[0]; /* TODO : Handle [a-] */
                if (regex[1] == '-' && regex[1] != ']') regex += 2;
                for (; i <= *regex; i++) {
                    if (b->s) b->s[b->l / BITS * ALPHABET + i] ^= 1 << (b->l % BITS);
                }
            } while (*++regex != ']');
            for (i = 0; inverse && b->s && i < ALPHABET; i++) {
                b->s[b->l / BITS * ALPHABET + i] ^= 1 << (b->l % BITS);
            }
        }
        else if (*regex == '.') {
            for (i = 0; i < ALPHABET && b->s; i++) {
                b->s[b->l / BITS * ALPHABET + i] |= 1 << (b->l % BITS);
            }
        }
        else if (*regex == '\\') {
            regex++; /* Strip off '\' */
            lastLen++;
            if (0 && *regex == 'E') {
                continue;
            }
            if (b->s) {
                if (strchr ("bBdDdSwW<>", *regex)) {
                    for (i = 0; i < ALPHABET; i++) {
                        if ((tolower (*regex) == 'b' ? (!isalnum (i) && i != '_') :
                             tolower (*regex) == 'd' ? isdigit (i) :
                             tolower (*regex) == 's' ? isspace (i) : 
                             isalnum (i) || i == '_') == (tolower (*regex) == *regex &&
                                                          *regex != '>')) b->s[b->l / BITS * ALPHABET + i] |=
                            1 << (b->l % BITS);
                    }
                }
                else b->s[b->l / BITS * ALPHABET + (
                                                    *regex == 'a' ? '\a' : *regex == 'e' ? '\x1b' : *regex == 'f' ?
                                                    '\f' : *regex == 'n' ? '\n' : *regex == 'r' ? '\r' : *regex == 
                                                    't' ? '\t' : *regex)] |= 1 << (b->l % BITS);
            }
        }
        else if (b->s) b->s[b->l / BITS * ALPHABET + *regex] |= 1 << (b->l%BITS);
        if (b->s) b->x[b->l / BITS] |= b->approxMode << (b->l % BITS);
        lastL = b->l - 1;
        b->l++;
        regex++;
    }
    if (infixEnd >= 0) NewJump (b, infixEnd, b->l - 1, infixRight);
    /* Jump from end of left hand to end */
    return regex - start;
}

int NewBitap (bitapType *b, const char *regex)
{ /* Returns a 0 or greater when successful, otherwise -x-1 where x is the
   * position of the first error. */
    int result;
    b->s = NULL;
    b->n = 0;
    b->l = 1; /* We create a 'start' state */
    result = ParseRegex (b, regex, NULL); /* First pass just to get sizes */
    if (result < 0) return result; 
    if (regex[result] != '\0') return -result - 1; /* Check for unmatched ')' */
    b->s = (int*) calloc ((b->l + BITS - 1) / BITS * ALPHABET, sizeof (*b->s));
    b->x = (int*) calloc ((b->l + BITS - 1) / BITS, sizeof (*b->x));
    b->j = (struct bitapJumpType *) malloc (b->n * sizeof (*b->j));
    b->n = 0;
    b->l = 1;
    b->x[0] = 1;
    b->approxMode = 1;
    ParseRegex (b, regex, NULL); /* No errors would occur. */
    return 0;
}

static void FromToForward (bitapType *b, int *r)
{ /* Handles the regular expression unglies at search time */
    int i;
    for (i = 0; i < b->n; i++) {
        if (r[b->j[i].foffset] & b->j[i].fmask) {
            r[b->j[i].toffset] |= b->j[i].tmask;
            
            if (b->j[i].t == infixRight) r[b->j[i].foffset] ^= b->j[i].fmask;
            /* Prevent flow from end of left hand side to start of right */
        }
    }
}

static void FromToReverse (bitapType *b, int *r)
{ /* Handles the regular expression unglies at search time */
    int i;
    for (i = b->n - 1; i >= 0; i--) {
        if (r[b->j[i].toffset] & b->j[i].tmask) {
            r[b->j[i].foffset] |= b->j[i].fmask;
            
            if (b->j[i].t == infixLeft) r[b->j[i].toffset] ^= b->j[i].tmask;
            /* Prevent flow from end of left hand side to start of right */
        }
    }
}


//        &b = the query word that had stuff done to it
//        in = the subject word
//        inl = the length of the subject word
//        e = max error
//        ereturn = errors in this word
//        breturn = the matched word
const char *FindWithBitap (bitapType *b, const char *in, int inl, int e, int *ereturn, const char **bReturn)
{ /* in is the 'haystack'. 'b' is the needle. e is the maximum of errors. */
    int i, j, carry, tmp, words = (b->l + BITS - 1) / BITS;
    int *r = (int*) calloc (words * (e + 1), sizeof (*r)), beste = e + 1;
    
    // this is what we return
    const char *best = NULL;
    
    
    for (i = 0; i <= e; i++) {
        r[i * words] = 1; /* Enable the start state. */
        FromToForward (b, r + i * words);
    }
    for (;; in++) 
    { /* Invarient : r[i].bit[j] implies r[i + words].bit[j] */
        for (i = 1; i <= e; i++) 
        { /* Now check insertion possibilities */
            for (carry = 0, j = 0; j < words; j++) 
            {
                tmp = r[(i - 1) * words + j];
                r[i * words + j] |= ((tmp<<1) | carry) & b->x[j];
                carry = tmp < 0 ? 1 : 0;
            }
            FromToForward (b, r + words * i);
        }
        while (beste && r[beste * words - 1] >> ((b->l - 1) % BITS)) 
        {
            beste--;
            best = in;
        }
        
        // if the word size is one
        if (!inl--) break;
        for (i = e; i >= 0; i--) 
        {
            /* The process is as follows, except that steps 3 and 4 are somewhat */
            /* combined */
            /* 1. Shift the last level ('e') */
            /* 2. Mask last level ('e') i.e. bits in the last level can only stay */
            /*    on if the character match the regex */
            /* 3. Use the on bits in the second last level to turn bits on in the */
            /*    last level : substitutions and deletions */
            /* 4. Do step 1 for the second last level */
            /* 5. Do step 2 for the second last level */
            /* ... */
            for (carry = 0, j = 0; j < words; j++) 
            {
                tmp = r[i * words + j];
                
                if (i < e) r[(i + 1) * words + j] |= tmp | (tmp << 1) | carry;
                /* tmp <=> deletions. */
                /* tmp << 1 | carry <=> substitutions */
                
                r[i * words + j] = ((tmp << 1) | carry) & b->s[*in + ALPHABET * j];
                carry = tmp < 0 ? 1 : 0;
            }
            if (i < e) FromToForward (b, r + words * (i + 1));
        }
        for (carry = 0, j = 0; j < words; j++) 
        { /* exact match row */
            /* tmp = r[j]; */
            r[j] &= b->s[*in + ALPHABET * j];
            /* carry = tmp < 0 ? 1 : 0; */
        }
        FromToForward (b, r);
    }
    
    // if you can't shift to the right here then there is no match
    if (!(r[words * (e + 1) - 1] >> ((b->l-1) % BITS))) best = NULL;
    
    // shift to the right as long as 'e' has a value and decrement at the end of the loop
    else while (e && r[words * e - 1] >> ((b->l - 1) % BITS)) e--;
    
    if (ereturn) *ereturn = e;
    
    if (bReturn && best) 
    {
        e = beste; /* These may not match if the pattern does not end with .**/
        memset (r, 0, sizeof (*r) * words * (e + 1));
        for (i = 0; i <= e; i++) 
        {
            r[(i + 1) * words - 1] |= 1 << ((b->l-1) % BITS);
            FromToReverse (b, r + i * words);
        }
        for (*bReturn = best; ; ) 
        {
            for (i = 1; i <= e; i++) 
            { /* Check for insertion possibilities */
                for (carry = 0, j = words - 1; j >= 0; j--) 
                {
                    tmp = r[(i - 1) * words + j];
                    r[i * words + j] |= (((unsigned) tmp >> 1) | carry) & b->x[j];
                    carry = tmp << (BITS - 1);
                }
                FromToReverse (b, r + words * i);
            }
            if (r[e * words] & 1) break;
            (*bReturn)--;
            for (i = e; i >= 0; i--) 
            {
                /* Same as above, except the mask takes place before the shift, */
                /* so split up steps 3 and 4 */
                if (i < e) {
                    for (carry = 0, j = words - 1; j >= 0; j--) 
                    {
                        tmp = r[i * words + j];
                        r[(i + 1) * words + j] |= tmp | ((unsigned) tmp >> 1) | carry;
                        carry = tmp << (BITS - 1);
                    }
                }
                for (carry = 0, j = words - 1; j >= 0; j--) 
                {
                    tmp = r[i * words + j] & b->s[**bReturn + ALPHABET * j];
                    r[i * words + j] = ((unsigned) tmp >> 1) | carry;
                    carry = tmp << (BITS - 1);
                }
                if (i < e) FromToReverse (b, r + words * (i + 1));
            }
            FromToReverse (b, r);
        } /* While searching for the beginning */
    } /* If we must return the beginning of the match */
    free (r);
    return best;
}

void DeleteBitap (bitapType *b)
{
    free (b->s);
    free (b->x);
    free (b->j);
}
