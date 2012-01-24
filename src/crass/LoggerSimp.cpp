// File: LoggerSimp.cpp
// Original Author: Michael Imelfort
// --------------------------------------------------------------------
//
// OVERVIEW:
// Implementation of LoggerSimp methods.
// No good for multithreaded apps
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

// system includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>

// local includes
#include "LoggerSimp.h"

LoggerSimp* LoggerSimp::mInstance = NULL;

LoggerSimp* LoggerSimp::Inst(void) {
    if(mInstance == NULL){
        mInstance = new LoggerSimp();
    }
    mInstance->setFileOpen(false);
    std::ofstream * fh = mInstance->getFhandle();
    fh = NULL;
    std::streambuf * buff = mInstance->getBuff();
    buff = NULL;
    return mInstance;
}

LoggerSimp::LoggerSimp() {}

LoggerSimp::~LoggerSimp(){
    if(mFileOpen)
        mInstance->closeLogFile();
    if(mTmpFH != NULL)
    {
    	delete mTmpFH;
    }
    if(mGlobalHandle != NULL)
    {
    	delete mGlobalHandle;
    }
    if(mInstance != NULL)
    {
        delete mInstance;
    }
}

void LoggerSimp::init(std::string logFile, int logLevel)
{
    mInstance->setLogLevel(logLevel);
    mInstance->setStartTime();
    if(logFile == "")
    {
        // set the logger to cout
        std::streambuf * buff = mInstance->getBuff();
        buff = std::cout.rdbuf();
        mGlobalHandle = new std::iostream(buff);
    }
    else
    {
        mInstance->setLogFile(logFile);
        mInstance->clearLogFile();
        mInstance->openLogFile();
    }
}

// Get methods
int LoggerSimp::getLogLevel(void)
{
    //-----
    // get the log level
    //
    return mLogLevel;
}

std::string LoggerSimp::getLogFile(void)
{
    //-----
    // the file we're logging to
    //
    return mLogFile;
}

bool LoggerSimp::isFileOpen(void)
{
    //-----
    // is the log file open?
    //
    return mFileOpen;
}

std::ofstream * LoggerSimp::getFhandle(void)
{
    //-----
    // get the fileHandle
    //
    return mFileHandle;
}

std::streambuf * LoggerSimp::getBuff(void)
{
    //-----
    // get the rbuff
    //
    return mBuff;
}

// Set Methods

void LoggerSimp::setFileOpen(bool isOpen)
{
    //----
    // set the file open flag
    //
    mFileOpen = isOpen;
}

void LoggerSimp::setStartTime(void)
{
    //----
    // set the start time
    //
    time ( &mStartTime );
}

void LoggerSimp::setLogLevel(int ll)
{
    //-----
    // set the log level
    //
    mLogLevel = ll;
}

void LoggerSimp::setLogFile(std::string lf)
{
    //-----
    // set the file name for the log file
    //
    mLogFile = lf;
}

// Operations
std::string LoggerSimp::int2Str(int input)
{
    //-----
    // curse c++ and their non toString()
    //
    std::stringstream ss;
    std::string s;
    ss << input;
    ss >> s;
    return s;
}

std::string LoggerSimp::timeToString(bool elapsed)
{
    //-----
    // get the time in a pretty form. Also can get time elapsed
    //
    struct tm * timeinfo;
    char buffer [80];
    
    time ( &mCurrentTime );
    
    if(elapsed)
    {
        std::string tmp = "";
        int tot_secs = (int)(difftime(mCurrentTime, mStartTime));
        int tot_days = tot_secs / 86400;
        if(tot_days)
        {
            tmp += int2Str(tot_days)+"d ";
            tot_secs = tot_secs - (tot_days * 86400);
        }
        int tot_hours = tot_secs / 3600;
        if(tot_hours)
        {
            tmp += int2Str(tot_hours)+"h ";
            tot_secs = tot_secs - (tot_hours * 3600);
        }
        int tot_mins = tot_secs / 60;
        if(tot_mins)
        {
            tmp += int2Str(tot_mins)+"m ";
            tot_secs = tot_secs - (tot_mins * 60);
        }
        tmp += int2Str(tot_secs)+"s";
        return tmp;
    }
    else
    {
        timeinfo = localtime ( &mCurrentTime );
        strftime (buffer,80,"%d/%m/%Y_%I:%M",timeinfo);
        std::string tmp(buffer);
        return tmp;
    }
}

void LoggerSimp::openLogFile(void)
{
    //-----
    // opens the log file
    //
	mTmpFH = mInstance->getFhandle();
	mTmpFH = new std::ofstream(getLogFile().c_str(), std::ios::app);
    mInstance->setFileOpen(true);
    std::streambuf * buff = mInstance->getBuff();
    buff = mTmpFH->rdbuf();
    mGlobalHandle = new std::iostream(buff);
}

void LoggerSimp::closeLogFile(void)
{
    //-----
    // closes the log file
    //
    std::ofstream * fh = mInstance->getFhandle();
    if(fh != NULL)
    {
        fh->close();
        delete fh;
    }
}

void LoggerSimp::clearLogFile(void)
{
    //-----
    // clears the log file
    //
    std::ofstream tmp_file(mLogFile.c_str(), std::ios::out);
    tmp_file.close();
}