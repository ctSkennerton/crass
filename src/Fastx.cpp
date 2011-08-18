/*
 *  Fastx.cpp is part of the crass project
 *  
 *  Created by Connor Skennerton on 10/08/11.
 *  Copyright 2011 Connor Skennerton. All rights reserved. 
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

#include <iostream>
#include <sstream>
#include <string>
#include "Fastx.h"


float Fastx::GCContent(void)
{
    std::string::iterator seq_iter = mSequence.begin();
    float count = 0.0;
    while (seq_iter != mSequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'G':
            case 'g':
            case 'C':
            case 'c':
                count++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return count / (float)mLength;
}


float Fastx::GCSkew(void)
{
    float g = 0.0;
    float c = 0.0;
    std::string::iterator seq_iter = mSequence.begin();
    while (seq_iter != mSequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'G':
            case 'g':
                g++;
                break;
            case 'C':
            case 'c':
                c++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return ((g-c)/(g+c));
}

float Fastx::ATContent(void)
{
    std::string::iterator seq_iter = mSequence.begin();
    float count = 0.0;
    while (seq_iter != mSequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'A':
            case 'a':
            case 'T':
            case 't':
                count++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return count / (float)mLength;
}


float Fastx::ATSkew(void)
{
    float a = 0.0;
    float t = 0.0;
    std::string::iterator seq_iter = mSequence.begin();
    while (seq_iter != mSequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'A':
            case 'a':
                a++;
                break;
            case 'T':
            case 't':
                t++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return ((a-t)/(a+t));
}

std::istream& Fasta::read(std::istream& stream)
{
    
    if (!(stream.good()))
    {
        stream.setstate (std::ios::badbit);
        return (stream);
    }
    
    if ( '>' != stream.get()) throw "badFormat, file not in FASTA format";
    
    mHeader.clear();
    mHeader.reserve(200);
    stream >> mHeader;
    mSequence.clear();
    mSequence.reserve(1000);
    stream.get();
    // check that there is some sequence associated with this header
    if(('>' == stream.peek()) || (stream.eof()))
    {
        throw "badFormat, no sequence in record";
    }
    
    // load the sequence into a stringstream
    std::stringstream ss;
    std::string tmp;
    do 
    {
        stream >> tmp;
        ss << tmp;
        stream.get();
    } while ('>' != stream.peek() && (!(stream.eof())));
    
    // copy the stringstream onto our sequence
    mSequence = ss.str();
    mLength = mSequence.length();
    return stream;
}

std::ostream& Fasta::print(std::ostream& s)
{
    s<<'>'<<mHeader<<std::endl<<mSequence;
    return s;
}


Fasta Fasta::subseq(int b, int l)
{
    Fasta tmp;
    
    tmp.mSequence = this->mSequence.substr(b,l);
    tmp.mHeader = this->mHeader;
    tmp.mLength = tmp.mSequence.length();
    
    return tmp;
    
}

Fasta Fasta::subseq(int b)
{
    Fasta tmp;
    tmp.mSequence = this->mSequence.substr(b);
    tmp.mHeader = this->mHeader;
    tmp.mLength = tmp.mSequence.length();
    return tmp;
}

Fasta Fasta::truncate(int b)
{
    Fasta tmp;
    tmp.mSequence = this->mSequence.substr(b);
    tmp.mHeader = this->mHeader;
    tmp.mLength = tmp.mSequence.length();
    return tmp;
}

Fasta Fasta::reverseComplement()
{
    Fasta tmp;
    tmp.header(this->header());
    std::stringstream revcomp_string;
    std::string::reverse_iterator rit;
    for ( rit = mSequence.rbegin() ; rit < mSequence.rend(); rit++ )
    {
		switch ((*rit)) 
        {
            case 'A':
            case 'a':
                revcomp_string << 'T';
                break;
            case 'C':
            case 'c':
                revcomp_string << 'G';
                break;
            case 'G':
            case 'g':
                revcomp_string << 'C';
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                revcomp_string << 'A';
                break;
            case 'N':
            case 'n':
                revcomp_string << 'N';
                break;
            case 'M':
            case 'm':
                revcomp_string <<'K'; 
                break;
            case 'R':
            case 'r':
                revcomp_string << 'Y';
                break;
            case 'W':
            case 'w':
                revcomp_string << 'W';
                break;
            case 'S':
            case 's':
                revcomp_string << 'S';
                break;
            case 'Y':
            case 'y':
                revcomp_string << 'R';
                break;
            case 'K':
            case 'k':
                revcomp_string << 'M';
                break;
            case 'V':
            case 'v':
                revcomp_string << 'B';
                break;
            case 'H':
            case 'h':
                revcomp_string << 'D';
                break;
            case 'D':
            case 'd':
                revcomp_string << 'H';
                break;
            case 'B':
            case 'b':
                revcomp_string << 'V';
                break;
            default:
                revcomp_string << 'N';
                break;
		}
	}
    tmp.seq(revcomp_string.str());
    return tmp;
}

Fasta Fasta::laurenize()
{
    Fasta seq2 = this->reverseComplement();
    if ((*this) < seq2)
    {
        return (*this);
    }
    return seq2;
}

std::istream& Fastq::read(std::istream& stream)
{
    
    if (!(stream.good()))
    {
        stream.setstate (std::ios::badbit);
        return (stream);
    }
    
    if ( '@' != stream.get()) throw "badFormat, file not in FASTQ format";
    
    mHeader.clear();
    mHeader.reserve(200);
    stream >> mHeader;
    
    mSequence.clear();
    mSequence.reserve(1000);
    
    // remove the trailing newline
    stream.get();
    // check that there is some sequence associated with this header
    if(('@' == stream.peek()) || ('+' == stream.peek()) || (stream.eof()))
    {
        throw "badFormat, no sequence in record";
    }
    
    // read in all of the lines of this sequence
    std::stringstream ss;
    std::string tmp;
    do 
    {
        stream >> tmp;
        ss << tmp;
        stream.get();
    } while (('+' != stream.peek()) && ('@' != stream.peek()) && (!(stream.eof())));
    
    mSequence = ss.str();
    
    //stream.get();
    if (('+' != stream.get()) || (stream.eof())) 
    {
        throw "badFormat, this record has no comment line";
    }
    
    mComment.clear();
    mComment.reserve(200);
    stream >> mComment;
    
    stream.get();
    // check that there is some quality associated with this header
    if(('@' == stream.peek()) || ('+' == stream.peek()) || (stream.eof()))
    {
        throw "badFormat, no quality in record";
    }
    
    mQuality.clear();
    mQuality.reserve(1000);
    // read in all of the lines of Quality scores for this sequence
    std::stringstream jj;
    std::string tp;
    do 
    {
        stream >> tp;
        jj << tp;
        stream.get();
    } while (('+' != stream.peek()) && ('@' != stream.peek()) && (!(stream.eof())));
    
    mQuality = jj.str();
    
    mLength = mSequence.length();
    return stream;
}

std::ostream& Fastq::print(std::ostream& s)
{
    s<<'@'<<mHeader<<std::endl<<mSequence<<std::endl<<'+'<<mComment<<std::endl<<mQuality;
    return s;
}

Fastq Fastq::subseq(int b, int l)
{
    Fastq tmp;
    
    tmp.mSequence = this->mSequence.substr(b,l);
    tmp.mHeader = this->mHeader;
    tmp.mComment = this->mComment;
    tmp.mQuality = this->mQuality.substr(b,l);
    tmp.mLength = tmp.mSequence.length();
    
    return tmp;
    
}

Fastq Fastq::subseq(int b)
{
    Fastq tmp;
    tmp.mSequence = this->mSequence.substr(b);
    tmp.mHeader = this->mHeader;
    tmp.mComment = this->mComment;
    tmp.mQuality = this->mQuality.substr(b);
    tmp.mLength = tmp.mSequence.length();
    return tmp;
}

Fastq Fastq::truncate(int b)
{
    Fastq tmp;
    tmp.mSequence = this->mSequence.substr(b);
    tmp.mHeader = this->mHeader;
    tmp.mComment = this->mComment;
    tmp.mQuality = this->mQuality.substr(b);
    tmp.mLength = tmp.mSequence.length();
    return tmp;
}

Fastq Fastq::reverseComplement()
{
    Fastq tmp;
    tmp.header(this->header());
    tmp.comment(this->comment());
    
    std::stringstream revcomp_string;
    std::stringstream rev_qaul;
    std::string::reverse_iterator rit = mSequence.rbegin();
    std::string::reverse_iterator qual_rit = mQuality.rbegin();
    
    while (rit != mSequence.rend()) 
    {
		switch ((*rit)) 
        {
            case 'A':
            case 'a':
                revcomp_string << 'T';
                break;
            case 'C':
            case 'c':
                revcomp_string << 'G';
                break;
            case 'G':
            case 'g':
                revcomp_string << 'C';
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                revcomp_string << 'A';
                break;
            case 'N':
            case 'n':
                revcomp_string << 'N';
                break;
            case 'M':
            case 'm':
                revcomp_string <<'K'; 
                break;
            case 'R':
            case 'r':
                revcomp_string << 'Y';
                break;
            case 'W':
            case 'w':
                revcomp_string << 'W';
                break;
            case 'S':
            case 's':
                revcomp_string << 'S';
                break;
            case 'Y':
            case 'y':
                revcomp_string << 'R';
                break;
            case 'K':
            case 'k':
                revcomp_string << 'M';
                break;
            case 'V':
            case 'v':
                revcomp_string << 'B';
                break;
            case 'H':
            case 'h':
                revcomp_string << 'D';
                break;
            case 'D':
            case 'd':
                revcomp_string << 'H';
                break;
            case 'B':
            case 'b':
                revcomp_string << 'V';
                break;
            default:
                revcomp_string << 'N';
                break;
		}
        rev_qaul << *qual_rit;
        rit++;
        qual_rit++;
	}
    tmp.seq(revcomp_string.str());
    tmp.qual(rev_qaul.str());
    return tmp;
}

Fastq Fastq::laurenize()
{
    Fastq seq2 = this->reverseComplement();
    if ((*this) < seq2)
    {
        return (*this);
    }
    return seq2;
}

// change the single ascii values into their corespnding Phread quality score
void Fastq::convertQualityToPhreadScore(int qualType)
{
    std::stringstream phred;
    std::string::iterator iter = mQuality.begin();
    short Q;
    while (iter != mQuality.end()) 
    {
        switch (qualType) 
        {
            case sanger:
                Q = (*iter) - 33;
                phred << Q;
                break;
            case illumina:
                Q = (*iter) - 64;
                phred << Q;
                break;
            default:
                break;
        }
        iter++;
    }
    mQuality = phred.str();

}
inline void Fastq::convertPhreadScoreToQuality(void)
{
    //convertQualityToPhreadScore(sanger);
}

void Fastq::convertPhreadScoreToQuality(int qualType)
{
    //std::stringstream qual;
    std::string::iterator iter = mQuality.begin();
    while (iter != mQuality.end()) 
    {
        //$q = chr(($Q<=93? $Q : 93) + 33);
        switch (qualType) 
        {
            case sanger:
                *iter = (char)(((*iter) <= 93 ? (*iter) : 93) + 33);
                break;
            case illumina:
                *iter = (char)(((*iter) <= 62 ? (*iter) : 62) + 64);
                break;
            default:
                break;
        }
        iter++;
    }

}

// change the single ascii values into their corespnding Phread quality score
// assuming that the quality scores are in sanger format
inline void Fastq::convertQualityToPhreadScore(void)
{
    convertQualityToPhreadScore(sanger);
}

// change the current ascii quality character into the corresponding Sanger type quality character
void Fastq::convertQualityToSangerAsciiValues(int qualType)
{
    std::string::iterator qual_iter;
    switch (qualType) 
    {
        case illumina:
            qual_iter = mQuality.begin();
            while (qual_iter != mQuality.end()) 
            {
                *qual_iter = *qual_iter - 31;
                qual_iter++;
            }
            break;
        case sanger:
            break;
        default:
            break;
    }

}

// change the current ascii quality character into the corresponding Sanger type quality character
// assumes that the current quality is in illumina format
inline void Fastq::convertQualityToSangerAsciiValues(void)
{
    convertQualityToSangerAsciiValues(illumina);
} 

// change the current quality values into their corresponding illumina values
void Fastq::convertQualityToIlluminaAsciiValues(int qualType)
{
    std::string::iterator qual_iter;
    switch (qualType) 
    {
        case sanger:
            qual_iter = mQuality.begin();
            while (qual_iter != mQuality.end()) 
            {
                *qual_iter = *qual_iter + 31;
                qual_iter++;
            }
            break;
        case illumina:
            break;
        default:
            break;
    }
}

inline void Fastq::convertQualityToIlluminaAsciiValues(void)
{
    convertQualityToIlluminaAsciiValues(sanger);
}

std::ostream & operator<< (std::ostream & s, Fastx & c)
{
    return c.print (s);
    
}

std::istream & operator>> (std::istream & s, Fastx &c)
{
    return c.read (s);
}

bool operator<(Fastx& f1, Fastx& f2)
{
    if (f1.mSequence < f2.mSequence) 
    {
        return true;
    }
    return false;
}

bool operator>(Fastx& f1, Fastx& f2)
{
    if (f1.mSequence > f2.mSequence) 
    {
        return true;
    }
    return false;
}