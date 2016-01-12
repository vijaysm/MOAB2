
MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest
AM_CPPFLAGS =  $(BGDEF)-DSRCDIR=$(srcdir) \
               $(BGDEF)-DMESHDIR=$(MESHDIR) \
               -I.. -I$(srcdir)/.. \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
	       -I$(top_builddir)/src/parallel 
            
TESTS = h5test \
        h5legacy \
        h5varlen \
        h5sets_test \
        h5regression \
        h5partial \
        h5portable

check_PROGRAMS = $(TESTS) dump_sets
LDADD = $(top_builddir)/src/libMOAB.la

h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
h5varlen_SOURCES = h5varlen.cpp
h5sets_test_SOURCES = h5sets_test.cpp
h5regression_SOURCES = h5regression.cpp
h5partial_SOURCES = h5partial.cpp
h5portable_SOURCES = h5portable.cpp

dump_sets_SOURCES = dump_sets.c
dump_sets_LDADD = $(top_builddir)/src/io/mhdf/libmhdf.la
