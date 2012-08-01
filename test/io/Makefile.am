LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/parallel

if NETCDF_FILE
  EXODUS_TEST = exodus_test
  NC_TEST = read_nc
else
  EXODUS_TEST = 
  NC_TEST = 
endif

if CCMIO_FILE
  CCMIO_TEST = ccmio_test
else
  CCMIO_TEST = 
endif

if HAVE_CGM
  CGM_TEST = read_cgm_test
else
  CGM_TEST = 
endif

TESTS = readutil_test \
	tqdcfr \
        vtk_test \
        cub_file_test \
	$(EXODUS_TEST) \
        gmsh_test \
	nastran_test \
	$(NC_TEST) \
        ideas_test \
        stl_test \
        smf_test \
	$(CCMIO_TEST) \
	$(CGM_TEST)

check_PROGRAMS = $(TESTS)

tqdcfr_SOURCES = $(top_srcdir)/src/io/Tqdcfr.cpp
tqdcfr_CXXFLAGS = -DTEST_TQDCFR $(CXXFLAGS)

#test_rms_SOURCES = test_rms.cpp
#test_exo_SOURCES = test_exo.cpp
vtk_test_SOURCES = $(srcdir)/../TestUtil.hpp VtkTest.cpp
cub_file_test_SOURCES = $(srcdir)/../TestUtil.hpp cub_file_test.cc
exodus_test_SOURCES = $(srcdir)/../TestUtil.hpp exodus_test.cc
gmsh_test_SOURCES = $(srcdir)/../TestUtil.hpp gmsh_test.cc
nastran_test_SOURCES = $(srcdir)/../TestUtil.hpp nastran_test.cc
read_nc_SOURCES = $(srcdir)/../TestUtil.hpp read_nc.cpp
ideas_test_SOURCES = $(srcdir)/../TestUtil.hpp ideas_test.cc
stl_test_SOURCES = $(srcdir)/../TestUtil.hpp stl_test.cc
smf_test_SOURCES = $(srcdir)/../TestUtil.hpp smf_test.cc
ccmio_test_SOURCES = $(srcdir)/../TestUtil.hpp ccmio_test.cpp
read_cgm_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_test.cpp
read_cgm_test_CXXFLAGS = $(CGM_CPPFLAGS)  $(CXXFLAGS)
readutil_test_SOURCES = $(srcdir)/../TestUtil.hpp readutil_test.cpp

MOSTLYCLEANFILES = dumped_acis.sat  \
                   mb_write.g \
                   mb_write2.g \
	           singlecyl_tmp.ccmg

