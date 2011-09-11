AM_CXXFLAGS = -Werror -pedantic -Wall -O2 -ggdb

AM_CPPFLAGS = 
if DEBUG
DEBUGFLAGS = -DDEBUG
AM_CPPFLAGS += $(DEBUGFLAGS)
endif

if SLOGGING
SLOGGINGFLAGS = -DSUPER_LOGGING
AM_CPPFLAGS += $(SLOGGINGFLAGS)
endif

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = crass

crass_SOURCES = bm.cpp Levensthein.cpp crass.cpp WuManber.cpp LoggerSimp.cpp SeqUtils.cpp CrisprNode.cpp NodeManager.cpp libcrispr.cpp WorkHorse.cpp SpacerInstance.cpp ReadHolder.cpp SmithWaterman.cpp StringCheck.cpp 

