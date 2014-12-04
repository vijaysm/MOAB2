MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src 

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/tools/mbcoupler \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/test 

lib_LTLIBRARIES = libmbcslam.la
LDADD = libmbcslam.la  $(top_builddir)/src/libMOAB.la

libmbcslam_la_SOURCES = \
   Intx2Mesh.cpp Intx2Mesh.hpp Intx2MeshOnSphere.cpp Intx2MeshOnSphere.hpp \
     CslamUtils.cpp CslamUtils.hpp Intx2MeshInPlane.cpp Intx2MeshInPlane.hpp IntxRllCssphere.hpp IntxRllCssphere.cpp
     
if PARALLEL     
    libmbcslam_la_SOURCES += wrap_intx.cpp 
endif

# these will be installed
libmbcslam_la_include_HEADERS = 

libmbcslam_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = intx_on_sphere_test  intx_in_plane_test spherical_area_test \
         case1_test  arc_intx_tests  IntxRllCssphere_test diffusion linear_remap
if PARALLEL
     TESTS += intx_mpas
endif

noinst_PROGRAMS =  proj1 create_dp 
if PARALLEL
  noinst_PROGRAMS += cslam_par_test intx_imesh
if ENABLE_FORTRAN
  noinst_PROGRAMS += advection
endif
endif

check_PROGRAMS = $(TESTS) 

intx_on_sphere_test_SOURCES = intx_on_sphere_test.cpp
intx_in_plane_test_SOURCES = intx_in_plane_test.cpp
# spec_visu_test_SOURCES = spec_visu_test.cpp
spherical_area_test_SOURCES = spherical_area_test.cpp
case1_test_SOURCES = case1_test.cpp
proj1_SOURCES = proj1.cpp
create_dp_SOURCES = create_dp.cpp
arc_intx_tests_SOURCES = arc_intx_tests.cpp
# this test is for netcdf only ; maybe convert to vtk for the test
IntxRllCssphere_test_SOURCES = IntxRllCssphere_test.cpp
linear_remap_SOURCES = linear_remap.cpp

if PARALLEL
  diffusion_SOURCES = diffusion.cpp
  intx_imesh_SOURCES = intx_imesh.cpp
  intx_imesh_LDADD = libmbcslam.la  $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la
  cslam_par_test_SOURCES = cslam_par_test.cpp
  intx_mpas_SOURCES = intx_mpas.cpp
if ENABLE_FORTRAN
  advection_SOURCES = advection.F90
  advection_LDADD = libmbcslam.la  $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la $(FCLIBS)
endif
endif 

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = intersect1.h5m \
                   lagr.h5m \
                   polyWithEdges.vtk \
                   lagr00.h5m \
                   intersect00.h5m \
                   lagrIni.h5m \
                   intx.vtk \
                   spectral.vtk \
                   intx1.vtk \
                   SpanEdges0.h5m \
                   SpanQuads0.h5m 
           
