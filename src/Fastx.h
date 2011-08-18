/*
 *  Fastx.h is part of the crass project
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

#ifndef crass_Fastx_h
#define crass_Fastx_h

#include <string>
class Fastx {
    
    
public:
    //members
    std::string mHeader;
    std::string mSequence;
    size_t      mLength;
    
    // constructor
    Fastx(void)
    {}
    
    // destructor
    virtual ~Fastx(void)
    {}
    
    // Getters & Setters
    inline std::string seq(void)
    {
        return mSequence;
    }
    inline void seq(std::string s)
    {
        mSequence = s;
    }
    
    inline void seq(const char * s)
    {
        mSequence = s;
    }    
    inline std::string header(void)
    {
        return mHeader;
    }
    
    inline void header(std::string h)
    {
        mHeader = h;
    }
 
    inline void header(const char * h)
    {
        mHeader = h;
    } 
    
    // member functions
    inline std::string substr(int begin, int length)
    {
        return mSequence.substr(begin, length);
    }
    
    inline std::string substr(int begin)
    {
        return mSequence.substr(begin);
    }
    
    inline char nucleotideAt(int ref)
    {
        return mSequence.at(ref);
    }
    
    inline char at(int ref)
    {
        return mSequence.at(ref);
    }
    
    float GCContent(void);
    
    inline int GCPercent(void)
    {
        return (int)(GCContent() * 100);
    }
    
    float GCSkew(void);
    
    float ATContent(void);
    
    inline int ATPercent(void)
    {
        return (int)(ATContent() * 100);
    }
    
    float ATSkew(void);
    
    // virtual functions
    virtual std::istream & read (std::istream & s) = 0;
    
    virtual std::ostream & print (std::ostream & s) = 0;
    
    virtual std::string qual(void) {return "";};
    
    virtual std::string comment(void) {return "";};
    
    virtual char qualityAt(int ref){return 0;};
    
    virtual void convertQualityToPhreadScore(void){};
    
    virtual void convertQualityToPhreadScore(int qualType){};
    
    virtual void convertPhreadScoreToQuality(void){};
    
    virtual void convertPhreadScoreToQuality(int qualType){};
    
    virtual void convertQualityToSangerAsciiValues(void){};
    
    virtual void convertQualityToSangerAsciiValues(int qualType){};
    
    virtual void convertQualityToIlluminaAsciiValues(void){};
    
    virtual void convertQualityToIlluminaAsciiValues(int qualType){};
    
    
    
};

class Fastq : public Fastx {
    
    
public:
    enum qualType {sanger, illumina};
    
    std::string mComment;
    std::string mQuality;
    
    Fastq(void)
    {
        mHeader = "";
        mSequence = "";
        mComment = "";
        mQuality = "";
        mLength = 0;
        
    }
    Fastq(std::string& h, std::string& s, std::string& c, std::string& q)
    {
        mHeader = h;
        mSequence = s;
        mComment = c;
        mQuality = q;
        mLength = s.length();
    }
    
    Fastq(const char * h, const char * s, const char * c, const char * q)
    {
        mHeader = h;
        mSequence = s;
        mComment = c;
        mQuality = q;
        mLength = strlen(s);
    }
    
    inline std::string qual(void)
    {
        return mQuality;
    }
    inline std::string comment(void)
    {
        return mComment;
    }
    inline void qual(std::string q)
    {
        mQuality = q;
    }
    inline void comment(std::string c)
    {
        mComment = c;
    }
    inline void qual(const char * q)
    {
        mQuality = q;
    }
    inline void comment(const char * c)
    {
        mComment = c;
    }
    inline char qualityAt(int ref)
    {
        return mQuality.at(ref);
    }
    
    Fastq reverseComplement(void);
    
    Fastq laurenize(void);
    
    Fastq subseq(int begin, int length);
    
    Fastq subseq(int begin);
    
    Fastq truncate(int begin);
    
    //Fasta to_fasta(void);
    
    std::istream& read(std::istream& s);
    
    std::ostream& print(std::ostream& s);
    
    inline void convertQualityToPhreadScore(void);
    
    void convertQualityToPhreadScore(int qualType);
    
    inline void convertPhreadScoreToQuality(void);

    void convertPhreadScoreToQuality(int qualType);

    inline void convertQualityToSangerAsciiValues(void);
    
    void convertQualityToSangerAsciiValues(int qualType);
    
    inline void convertQualityToIlluminaAsciiValues(void);
    
    void convertQualityToIlluminaAsciiValues(int qualType);
    
    
};


class Fasta : public Fastx {
    
    
public:
    Fasta(void)
    {
        mHeader = "";
        mSequence = "";
        mLength = 0;
    }
    
    Fasta(std::string& h, std::string& s)
    {
        mHeader = h;
        mSequence = s;
        mLength = s.length();
    }
    Fasta(const char * h, const char * s)
    {
        mHeader = h;
        mSequence = s;
        mLength = strlen(s);
    }  
    
    ~Fasta(void)
    {
        mHeader.clear();
        mSequence.clear();
    }
    
    Fasta reverseComplement(void);
    
    Fasta laurenize(void);
    
    Fasta subseq(int begin, int length);
    
    Fasta subseq(int begin);
    
    Fasta truncate(int begin);
    
    std::istream& read(std::istream& s);
    
    std::ostream& print(std::ostream& s);
    
    
    
    
};

std::ostream& operator<< (std::ostream& s, Fastx& c);

std::istream& operator>> (std::istream& s, Fastx& c);

bool operator<(Fastx& f1, Fastx& f2);

bool operator>(Fastx& f1, Fastx& f2);


#endif
