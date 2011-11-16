AM_CXXFLAGS = -Werror -pedantic -Wall -g

AM_LDFLAGS = @XERCES_LDFLAGS@
if ASSEMBLY_WRAPPER
    LIBS += -lxerces-c
endif
if DEBUG
    AM_CXXFLAGS += -ggdb
endif
if MEMCHECK
    AM_CXXFLAGS +=  -O0
else
    AM_CXXFLAGS += -O2
endif


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = crass

crass_SOURCES =\
PatternMatcher.cpp PatternMatcher.h\
Rainbow.cpp Rainbow.h\
crass.cpp crass.h\
WuManber.cpp WuManber.h\
LoggerSimp.cpp LoggerSimp.h\
SeqUtils.cpp SeqUtils.h\
CrisprNode.cpp CrisprNode.h\
NodeManager.cpp NodeManager.h\
libcrispr.cpp libcrispr.h\
WorkHorse.cpp WorkHorse.h\
SpacerInstance.cpp SpacerInstance.h\
ReadHolder.cpp ReadHolder.h\
SmithWaterman.cpp SmithWaterman.h\
StringCheck.cpp StringCheck.h\
kseq.h\
GraphDrawingDefines.h\
crassDefines.h\
StlExt.h

if ASSEMBLY_WRAPPER
crass_SOURCES +=\
CrassXML.cpp CrassXML.h\
AssemblyWrapper.cpp AssemblyWrapper.h
endif

