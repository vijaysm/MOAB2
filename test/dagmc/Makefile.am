LDADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/tools/dagmc/DagMC.lo

MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/tools/dagmc 

TESTS = dagmc_test

check_PROGRAMS = $(TESTS)

dagmc_test_SOURCES = $(srcdir)/../TestUtil.hpp dagmc_test.cpp
dagmc_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CXXFLAGS) -lcgm -lcubiti19


