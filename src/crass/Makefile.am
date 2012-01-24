AM_CXXFLAGS = -Werror -pedantic -Wall


if DEBUG
    AM_CXXFLAGS += -ggdb -O2
endif
if MEMCHECK
    AM_CXXFLAGS +=  -g -O0
endif


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = crass
INCLUDES = -I$(top_srcdir)/src/libcrispr
crass_LDADD = $(top_builddir)/src/libcrispr/libcrispr.la
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
StlExt.h\
StatsManager.h

if ASSEMBLY_WRAPPER
crass_SOURCES +=\
AssemblyWrapper.cpp AssemblyWrapper.h
endif

