/*
 *  Search.cpp
 *  Created by Connor Skennerton.
 *  Copyright 2011-2013 Connor Skennerton & Michael Imelfort. All rights reserved.
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

#ifndef __crass__Search__
#define __crass__Search__

#include <cstdio>
#include <set>
#include <string>
#include "RepeatArray.h"
#include "Sequence.h"
// --------------------------------------------------------------------
// SEARCH ALGORITHM PARAMETERS
// --------------------------------------------------------------------
#define CRASS_DEF_MIN_SEARCH_WINDOW_LENGTH         (6)
#define CRASS_DEF_MAX_SEARCH_WINDOW_LENGTH         (9)
#define CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH     (8)
#define CRASS_DEF_SCAN_LENGTH                      (30)
#define CRASS_DEF_SCAN_CONFIDENCE                  (0.70)
#define CRASS_DEF_TRIM_EXTEND_CONFIDENCE           (0.5)
#define CRASS_DEF_MIN_DR_SIZE                   (23)                  // minimum direct repeat size
#define CRASS_DEF_MAX_DR_SIZE                   (47)                  // maximum direct repeat size
#define CRASS_DEF_MIN_SPACER_SIZE               (26)                  // minimum spacer size
#define CRASS_DEF_MAX_SPACER_SIZE               (50)                  // maximum spacer size
// --------------------------------------------------------------------
// HARD CODED PARAMS FOR DR FILTERING
// --------------------------------------------------------------------
#define CRASS_DEF_LOW_COMPLEXITY_THRESHHOLD        (0.75)
#define CRASS_DEF_SPACER_OR_REPEAT_MAX_SIMILARITY  (0.82)
#define CRASS_DEF_SPACER_TO_SPACER_LENGTH_DIFF     (12)
#define CRASS_DEF_SPACER_TO_REPEAT_LENGTH_DIFF     (30)
#define CRASS_DEF_DEFAULT_MIN_NUM_REPEATS          (3)
#define CRASS_DEF_KMER_MAX_ABUNDANCE_CUTOFF	       (0.23)			// DRs should NOT have kmers more abundant than this percentage!

#define logInfo(cOUTsTRING, lEVEL) std::cerr <<__FILE__<<" "<<__LINE__<<">> "<< cOUTsTRING << std::endl;
namespace crass {
    
class Search {
  

public:
    Search() : 
            mMinRepeatLength(CRASS_DEF_MIN_DR_SIZE),
            mMaxRepeatLength(CRASS_DEF_MAX_DR_SIZE),
            mMinSpacerLength(CRASS_DEF_MIN_SPACER_SIZE),
            mMaxSpacerLength(CRASS_DEF_MAX_SPACER_SIZE),
            mSearchWindowLength(CRASS_DEF_OPTIMAL_SEARCH_WINDOW_LENGTH),
            mScanRange(CRASS_DEF_SCAN_LENGTH),
            mMinSeedCount(CRASS_DEF_DEFAULT_MIN_NUM_REPEATS)
    {}
    
    // get/set
    int minRepeatLength() {
        return mMinRepeatLength;
    }
    void minRepeatLength(int i) {
        mMinRepeatLength = i;
    }
    int maxRepeatLength() {
        return mMaxRepeatLength;
    }
    void maxRepeatLength(int i) {
        mMaxRepeatLength = i;
    }
    int minSpacerLength() {
        return mMinSpacerLength;
    }
    void minSpacerLength(int i) {
        mMinSpacerLength = i;
    }
    int maxSpacerLength() {
        return mMaxSpacerLength;
    }
    void maxSpacerLength(int i) {
        mMaxSpacerLength = i;
    }
#if 0
    // for single reads/files
    bool searchFileSerial(FILE *file);
    
    bool searchFileParallel(FILE *file);
    
    bool searchFileForPatterns(FILE *file, crass::list& patters);
    bool searchFileForPatterns(FILE *file, crass::list& patters, std::set<std::string>& skipReads);
    
    bool search(std::string& read, crass::RepeatArray& drPositions);
    
    //for read pairs/file pairs
    bool searchFileSerial(FILE *file1, FILE *file2);
    
    bool searchFileParallel(FILE *file, FILE *file2);
    
    bool searchFileForPatterns(FILE *file1, FILE *file2, crass::list& patters);
    bool searchFileForPatterns(FILE *file1, FILE *file2, crass::list& patters, std::set<std::string>& skipReads);
    
    bool search(std::string& read1, std::string& read2, crass::RepeatArray& drPositions);
#endif
    int longReadSearch(crass::RawRead& read);
    
    int scanRight(crass::RawRead& read,
                  std::string& pattern);
    
    int extendPreRepeat(crass::RawRead& read);
    
private:
    bool qcFoundRepeats(crass::RawRead& read);
    
    bool isRepeatLowComplexity(std::string& repeat);
    
    bool doesRepeatHaveHighlyAbundantKmers(std::string& directRepeat, float& max_count);
    
    bool doesRepeatHaveHighlyAbundantKmers(std::string& directRepeat);

#if 0
    bool searchBrute(crass::RawRead& read);
    bool searchCRT(crass::RawRead& read);
    
    void addReadHolder(ReadHolder& tmp_holder);
    
    // members
    crass::Storage *mReadStore;
#endif
    int mMinRepeatLength;
    int mMaxRepeatLength;
    int mMinSpacerLength;
    int mMaxSpacerLength;
    int mSearchWindowLength;
    int mScanRange;
    int mMinSeedCount;
    
};
    
}
#endif /* defined(__crass__Search__) */
