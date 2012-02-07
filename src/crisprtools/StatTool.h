/*
 * StatTool.h
 *
 * Copyright (C) 2011 - Connor Skennerton
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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

#ifndef STATTOOL_H
#define STATTOOL_H 

#include <vector>
#include <string>
#include <set>
#include "Xml.h"
#include "StlExt.h"


#define SPACER_CHAR '+'
#define FLANKER_CHAR '~'
#define REPEAT_CHAR '-'
typedef struct __AStats {
    
    int total_groups;
    int total_spacers;
    int total_dr;
    int total_flanker;
    int total_spacer_length;
    int total_spacer_cov;
    int total_dr_length;
    int total_flanker_length;
    } AStats;

class StatManager {
    std::vector<int> SM_SpacerLength;
    std::vector<int> SM_SpacerCoverage;
    std::vector<int> SM_RepeatLength;
    std::vector<int> SM_FlankerLength;
    int SM_SpacerCount;
    int SM_RepeatCount;
    int SM_FlankerCount;
    std::string SM_ConsensusRepeat;
    std::string SM_Gid;
    
public:
    
    StatManager() {
        SM_FlankerCount = 0;
        SM_RepeatCount = 0;
        SM_SpacerCount = 0;
    }
    
    inline std::string getConcensus(void){return SM_ConsensusRepeat;}
    inline std::string getGid(void){return SM_Gid;}
    
    inline int getSpacerCount(void) {return SM_SpacerCount;}
    inline int getRpeatCount(void) {return SM_RepeatCount;}
    inline int getFlankerCount(void) {return SM_FlankerCount;}
    

    
    inline void incrementSpacerCount(void) {++SM_SpacerCount;}
    inline void incrementRpeatCount(void) {++SM_RepeatCount;}
    inline void incrementFlankerCount(void) {++SM_FlankerCount;}
    
    inline std::vector<int> getSpLenVec(void){return SM_SpacerLength;}
    inline std::vector<int> getSpCovVec(void){return SM_SpacerCoverage;}
    inline std::vector<int> getRepLenVec(void){return SM_RepeatLength;}
    inline std::vector<int> getFlLenVec(void){return SM_FlankerLength;}
    
    
    
    inline void setConcensus(std::string s){ SM_ConsensusRepeat = s;}
    inline void setGid(std::string g){ SM_Gid = g;}
    
    inline void setSpacerCount(int i) { SM_SpacerCount = i;}
    inline void setRepeatCount(int i) { SM_RepeatCount = i;}
    inline void setFlankerCount(int i) { SM_FlankerCount = i;}
    
    // mode
    
    inline int modeSpacerL(void){return mode(SM_SpacerLength);}
    inline int modeSpacerC(void){return mode(SM_SpacerCoverage);}
    inline int modeRepeatL(void){return mode(SM_RepeatLength);}
    inline int modeFlankerL(void){return mode(SM_FlankerLength);}
    
    // median;
    
    inline int medianSpacerL(void){return median(SM_SpacerLength);}
    inline int medianSpacerC(void){return median(SM_SpacerCoverage);}
    inline int medianRepeatL(void){return median(SM_RepeatLength);}
    inline int medianFlankerL(void){return median(SM_FlankerLength);}
    
    // mean
    
    inline int meanSpacerL(void) {return mean(SM_SpacerLength);}
    inline int meanSpacerC(void) {return mean(SM_SpacerCoverage);}
    inline int meanRepeatL(void) {return mean(SM_RepeatLength);}
    inline int meanFlankerL(void) {return mean(SM_FlankerLength);}
    
    inline void addSpLenVec(int i ){return SM_SpacerLength.push_back(i);}
    inline void addSpCovVec(int i ){return SM_SpacerCoverage.push_back(i);}
    inline void addRepLenVec(int i ){return SM_RepeatLength.push_back(i);}
    inline void addFlLenVec(int i ){return SM_FlankerLength.push_back(i);}
    
};

class StatTool {

    std::set<std::string> ST_Groups;
    
    std::vector<StatManager *> ST_StatsVec;
    
    bool ST_Pretty;
    bool ST_AssemblyStats;
    bool ST_Subset;
    std::string ST_OutputFileName;
    bool ST_WithHeader;
    bool ST_AggregateStats;
    bool ST_Tabular;
    std::string ST_Separator;
    
public:
    StatTool() {

        ST_Pretty = false;
        ST_Subset = false;
        ST_AssemblyStats = false;
        ST_WithHeader = false;
        ST_AggregateStats = false;
        ST_Tabular = true;
        ST_Separator = "\t";
    }
    ~StatTool();


    
    void generateGroupsFromString(std::string str);
    int processOptions(int argc, char ** argv);
    int processInputFile(const char * inputFile);
    void parseGroup(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void parseData(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager);
    void parseDrs(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager);
    void parseSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager);
    void parseFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser, StatManager * statManager);
    
//    void parseAssembly(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
//    void parseContig(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
//    void parseCSpacer(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
//    void parseLinkSpacers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
//    void parseLinkFlankers(xercesc::DOMElement * parentNode, crispr::XML& xmlParser);
    void calculateAgregateSTats(AStats * agregateStats);
    void prettyPrint(StatManager * sm);
    void printHeader(void);
    void printTabular(StatManager * sm);
    void printAggregate(AStats * agregateStats);
    std::vector<StatManager *>::iterator begin(){return ST_StatsVec.begin();}
    std::vector<StatManager *>::iterator end(){return ST_StatsVec.end();}


};
int statMain(int argc, char ** argv);
void statUsage(void);
#endif
