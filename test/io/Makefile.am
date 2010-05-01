LDADD = $(top_builddir)/src/libMOAB.la

EXTRA_DIST = brick_cubit10.1.cub \
             brick_cubit10.2.cub \
             brick_cubit10.cub \
	     gmsh2.msh \
             ho_test.g \
             ho_test.cub \
             test.cub \
             test.nas \
             test.unv \
	     sample.stl \
             three.smf


AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/parallel

if NETCDF_FILE
  EXODUS_TEST = exodus_test
else
  EXODUS_TEST = 
endif

TESTS = tqdcfr \
        vtk_test \
        cub_file_test \
	$(EXODUS_TEST) \
        gmsh_test \
	nastran_test \
        ideas_test \
        stl_test \
        smf_test

check_PROGRAMS = $(TESTS)

tqdcfr_SOURCES = $(top_srcdir)/src/io/Tqdcfr.cpp
tqdcfr_CXXFLAGS = -DTEST_TQDCFR $(CXXFLAGS)
tqdcfr_DEPENDENCIES = $(LDADD)

#test_rms_SOURCES = test_rms.cpp
#test_rms_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_exo_SOURCES = test_exo.cpp
#test_exo_DEPENDENCIES = $(top_builddir)/libMOAB.la

vtk_test_SOURCES = $(srcdir)/../TestUtil.hpp VtkTest.cpp
vtk_test_DEPENDENCIES = $(LDADD)

cub_file_test_SOURCES = $(srcdir)/../TestUtil.hpp cub_file_test.cc
cub_file_test_DEPENDENCIES = $(LDADD)

exodus_test_SOURCES = $(srcdir)/../TestUtil.hpp exodus_test.cc
exodus_test_DEPENDENCIES = $(LDADD)

gmsh_test_SOURCES = $(srcdir)/../TestUtil.hpp gmsh_test.cc
gmsh_test_DEPENDENCIES = $(LDADD)

nastran_test_SOURCES = $(srcdir)/../TestUtil.hpp nastran_test.cc
nastran_test_DEPENDENCIES = $(LDADD)

ideas_test_SOURCES = $(srcdir)/../TestUtil.hpp ideas_test.cc
ideas_test_DEPENDENCIES = $(LDADD)

stl_test_SOURCES = $(srcdir)/../TestUtil.hpp stl_test.cc
stl_test_DEPENDENCIES = $(LDADD)

smf_test_SOURCES = $(srcdir)/../TestUtil.hpp smf_test.cc
smf_test_DEPENDENCIES = $(LDADD)

MOSTLYCLEANFILES = dumped_acis.sat 

