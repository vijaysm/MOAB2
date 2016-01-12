LDADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/tools/dagmc/DagMC.lo

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/tools/dagmc 
TESTS = 
if HAVE_HDF5
 TESTS += dagmc_simple_test \
        dagmc_rayfire_test \
        dagmc_pointinvol_test
endif
check_PROGRAMS = $(TESTS)

dagmc_simple_test_SOURCES = $(srcdir)/../TestUtil.hpp dagmc_simple_test.cpp
dagmc_simple_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CXXFLAGS) $(CGM_LIBS)
dagmc_rayfire_test_SOURCES = $(srcdir)/../TestUtil.hpp dagmc_rayfire_test.cpp
dagmc_rayfire_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CXXFLAGS) $(CGM_LIBS)
dagmc_pointinvol_test_SOURCES = $(srcdir)/../TestUtil.hpp dagmc_pointinvol_test.cpp
dagmc_pointinvol_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CXXFLAGS) $(CGM_LIBS)


