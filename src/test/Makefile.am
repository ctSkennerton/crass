TESTS = crass-test
check_PROGRAMS = crass-test
AM_CXXFLAGS = -I$(top_builddir)/src/crass/
AM_LDFLAGS = @zlib_flags@
crass_test_SOURCES = \
test_readholder.cpp\
test_libcrispr.cpp\
test_main.cpp

crass_test_LDADD = $(top_builddir)/src/crass/libcrass.a $(top_builddir)/src/aho-corasick/libacism.a
