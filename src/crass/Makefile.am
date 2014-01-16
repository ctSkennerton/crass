# this lists the binaries to produce, the (non-PHONY, binary) targets in
bin_PROGRAMS = crass


AM_CXXFLAGS = -std=c++0x -pedantic -Wall -DNO_OPENMP
AM_LDFLAGS =  @zlib_flags@



crass_SOURCES =\
PatternMatcher.cpp\
kseq.cpp Search.cpp\
Sequence.cpp\
RepeatArray.cpp\
Storage.cpp\
main.cpp\
LoggerSimp.cpp\
StringCheck.cpp\
SeqUtils.cpp\
cdhit-common.cpp\
cdhit-est.cpp\
WuManber.cpp\
Graph.cpp
