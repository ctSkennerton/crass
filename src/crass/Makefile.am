


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = crass

INCLUDES = -I$(top_srcdir)/src/libcrispr
crass_LDADD = $(top_builddir)/src/libcrispr/libcrispr.la

crass_CXXFLAGS = @XERCES_CPPFLAGS@ -Werror -pedantic -Wall
crass_LDFLAGS = @zlib_flags@


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
StatsManager.h\
SearchChecker.cpp\
SearchChecker.h

if ASSEMBLY_WRAPPER
crass_SOURCES +=\
AssemblyWrapper.cpp AssemblyWrapper.h
endif

