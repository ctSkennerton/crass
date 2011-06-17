//            File: LoggerSimp.h
// Original Author: Michael Imelfort
// --------------------------------------------------------------------
//
// OVERVIEW:
// This file contains the class definition for a simple output logger
// no good for multi threaded apps
//
// This is for runtime logging. For compile time and paranoid logging see
// the file: paranoid.h
//
// --------------------------------------------------------------------
// Copyright (C) 2009 2010 2011 Michael Imelfort and Dominic Eales
//
// This file is part of the Sassy Assembler Project.
//
// Sassy is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Sassy is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Sassy.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
//
//                        A
//                       A B
//                      A B R
//                     A B R A
//                    A B R A C
//                   A B R A C A
//                  A B R A C A D
//                 A B R A C A D A
//                A B R A C A D A B 
//               A B R A C A D A B R  
//              A B R A C A D A B R A 
//

#ifndef LoggerSimp_h
#define LoggerSimp_h

#include <time.h>
#include <iostream>
#include "crass_defines.h"
using namespace std;

// for making the main logger
#define intialiseGlobalLogger(lOGfILE, lOGlEVEL) logger->init(lOGfILE, lOGlEVEL)

// for determining if logging is possible at a given level
#define willLog(lOGlEVEL) (logger->getLogLevel() >= lOGlEVEL)

class LoggerSimp {
public:
    
    // Constructor/Destructor
    static LoggerSimp* Inst(void);                                  // get the singleton logger
    ~LoggerSimp();                                                  // kill it! [call this at the end of main()]
    
    // Get methods
    int getLogLevel(void);                                          // get the log level
    std::string getLogFile(void);                                        // the file we're logging to
    bool isFileOpen(void);                                          // is the log file open?
    std::ofstream * getFhandle(void);                                    // get the fileHandle
    std::streambuf * getBuff(void);                                      // get the rbuff
    
    // Set Methods
    void init(std::string logFile, int logLevel);                        // don't call this directly, use the macro instead
    void setStartTime(void);                                        // what time was the logger instantiated?
    void setLogLevel(int ll);                                       // set the log level
    void setLogFile(std::string lf);                                     // set the file name for the log file
    void setFileOpen(bool isOpen);                                  // set if the file is open
    
    // Operations
    std::string int2Str(int input);                                      // convert an into to a string
    std::string timeToString(bool elapsed);                              // write out the current time, prettylike
    void closeLogFile(void);                                        // close the log file down
    void openLogFile(void);                                         // open the log file
    void clearLogFile(void);                                        // clear the logFile at the start
    
    std::iostream * mGlobalHandle;                                       // what we realy write to
    
protected:
    LoggerSimp();
    
private:
    
    static LoggerSimp * mInstance;                                  // the internal instance for the singleton
    
    std::ofstream * mFileHandle;                                         // for writing to files
    std::streambuf *mBuff;                                               // for holding rbuffs
    
    std::string mLogFile;                                                // this is the file we'll be writing out to
    int mLogLevel;                                                  // which logging level are we at?
    time_t mStartTime;                                              // the time when the logger was created
    time_t mCurrentTime;                                            // now, .. no ... NOW! NOW!
    bool mFileOpen;                                                 // is the log file open?
};

static LoggerSimp* logger = LoggerSimp::Inst();                     // this makes the singleton available to all classes
// which include LoggerSimp.h
// for logging info
#define logInfo(cOUTsTRING, ll) { \
if(logger->getLogLevel() >= ll) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tI   " << cOUTsTRING << std::endl; \
} \
}

// for errors
#define logError(cOUTsTRING) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tERR " << __FILE__ << " : " << __PRETTY_FUNCTION__ << " : " << __LINE__ << ": " <<  cOUTsTRING << std::endl; \
}

// for warnings
#define logWarn(cOUTsTRING, ll) { \
if(logger->getLogLevel() >= ll) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tW   " << cOUTsTRING << std::endl; \
} \
}

// time stamp
#define logTimeStamp() { \
(*(logger->mGlobalHandle)) << "----------------------------------------------------------------------\n----------------------------------------------------------------------\n -- " << logger->timeToString(false) << "  --  " << LONG_NAME<<" ("<<PRG_NAME<<")" << " --  Version: " << MCD_VERSION << " --\n----------------------------------------------------------------------\n----------------------------------------------------------------------\n" << std::endl; \
}

#ifdef SUPER_LOGGING
#undef logInfo
#undef logError
#undef logWarn

// for logging info
#define logInfo(cOUTsTRING, ll) { \
if(logger->getLogLevel() >= ll) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tI   " << __FILE__ << " : " << __PRETTY_FUNCTION__ << " : " << __LINE__ << ": " <<  cOUTsTRING << std::endl; \
} \
}

// for errors
#define logError(cOUTsTRING) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tERR " << __FILE__ << " : " << __PRETTY_FUNCTION__ << " : " << __LINE__ << ": " <<  cOUTsTRING << std::endl; \
}

// for warnings
#define logWarn(cOUTsTRING, ll) { \
if(logger->getLogLevel() >= ll) { \
(*(logger->mGlobalHandle)) << logger->timeToString(true) << "\tW   " << __FILE__ << " : " << __PRETTY_FUNCTION__ << " : " << __LINE__ << ": " <<  cOUTsTRING << std::endl; \
} \
}

#endif // __SUPER_LOGGING

#endif //LoggerSimp_h