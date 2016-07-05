# this lists the binaries to produce, the (non-PHONY, binary) targets in
bin_PROGRAMS = crass crisprtools
noinst_LIBRARIES = libcrass.a
if ASSEMBLY_WRAPPER
bin_PROGRAMS += crass-assembler
endif

AM_CXXFLAGS = @XERCES_CPPFLAGS@ -pedantic -Wall
AM_LDFLAGS = @XERCES_LDFLAGS@ @zlib_flags@ @XERCES_LIBS@

crass_LDADD = libcrass.a $(top_builddir)/src/aho-corasick/libacism.a

crisprtools_LDADD = @GV_LIBS@ 

libcrass_a_SOURCES =\
PatternMatcher.cpp PatternMatcher.h\
Rainbow.cpp Rainbow.h\
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
kseq.cpp kseq.h\
GraphDrawingDefines.h\
crassDefines.h\
StatsManager.h\
SearchChecker.cpp SearchChecker.h\
ksw.c ksw.h\
Types.h\
Aligner.cpp Aligner.h\
base.cpp\
parser.cpp\
reader.cpp\
writer.cpp\
 $(top_builddir)/config.h

crass_SOURCES =\
crass.cpp crass.h

crass_assembler_SOURCES =\
AssemblyWrapper.cpp AssemblyWrapper.h\
$(top_builddir)/config.h\
crassDefines.h\
kseq.cpp kseq.h\
SeqUtils.cpp SeqUtils.h\
base.cpp\
parser.cpp\
reader.cpp\
writer.cpp

crisprtools_SOURCES = \
	crisprtools.cpp \
	MergeTool.cpp \
	MergeTool.h \
	SplitTool.cpp \
	SplitTool.h \
	FilterTool.cpp \
	FilterTool.h \
	ExtractTool.cpp \
	ExtractTool.h \
	SanitiseTool.cpp \
	SanitiseTool.h \
	StatTool.cpp \
	StatTool.h \
	Utils.cpp \
	Utils.h \
	Rainbow.cpp \
	Rainbow.h \
	RemoveTool.h \
	RemoveTool.cpp \
base.cpp\
parser.cpp\
reader.cpp\
writer.cpp

if FOUND_GRAPHVIZ_LIBRARIES
crisprtools_SOURCES += DrawTool.cpp DrawTool.h CrisprGraph.cpp CrisprGraph.h
endif

