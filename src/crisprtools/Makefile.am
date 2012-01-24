## Process this file with automake to produce Makefile.in

INCLUDES = \
	-I$(top_srcdir)/src/libcrispr

bin_PROGRAMS = crisprtools

crisprtools_LDADD = $(top_builddir)/src/libcrispr/libcrispr.la
crisprtools_LDFLAGS = @GRAPHVIZ_LIBRARIES@
crisprtools_SOURCES = \
	main.cpp \
	StlExt.h \
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
	RemoveTool.cpp 
    
if FOUND_GRAPHVIZ_LIBRARIES
crisprtools_SOURCES += DrawTool.cpp DrawTool.h CrisprGraph.cpp CrisprGraph.h
endif



