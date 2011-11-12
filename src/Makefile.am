AM_CXXFLAGS = -Werror -pedantic -Wall 
LIBS += @LIBXERCES_C@
if DEBUG
    AM_CXXFLAGS += -ggdb
endif


# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = crass
 
if ASSEMBLY_WRAPPER
crass_SOURCES = PatternMatcher.cpp Rainbow.cpp crass.cpp WuManber.cpp LoggerSimp.cpp SeqUtils.cpp CrisprNode.cpp NodeManager.cpp libcrispr.cpp WorkHorse.cpp SpacerInstance.cpp ReadHolder.cpp SmithWaterman.cpp StringCheck.cpp CrassXML.cpp AssemblyWrapper.cpp
else
crass_SOURCES = PatternMatcher.cpp Rainbow.cpp crass.cpp WuManber.cpp LoggerSimp.cpp SeqUtils.cpp CrisprNode.cpp NodeManager.cpp libcrispr.cpp WorkHorse.cpp SpacerInstance.cpp ReadHolder.cpp SmithWaterman.cpp StringCheck.cpp
endif
