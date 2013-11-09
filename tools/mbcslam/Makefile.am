MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src 

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/moab/point_locater/lotte \
               -I$(top_srcdir)/tools/mbcoupler \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/test 
if ENABLE_mbzoltan
  AM_CPPFLAGS +=  $(ZOLTAN_INC_FLAGS) \
               -I$(top_srcdir)/tools/mbzoltan \
               -I$(top_builddir)/tools/mbzoltan
if HAVE_CGM
  AM_CPPFLAGS += @CGM_CPPFLAGS@ -DCGM @MOAB_CGM_DEFINES@
endif
endif

lib_LTLIBRARIES = libmbcslam.la
libmbcslam_la_LIBADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/itaps/imesh/libiMesh.la \
         $(top_builddir)/tools/mbcoupler/libmbcoupler.la
         
LDADD = $(top_builddir)/src/libMOAB.la libmbcslam.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la

libmbcslam_la_SOURCES = \
   Intx2Mesh.cpp Intx2Mesh.hpp Intx2MeshOnSphere.cpp Intx2MeshOnSphere.hpp \
     CslamUtils.cpp CslamUtils.hpp Intx2MeshInPlane.cpp Intx2MeshInPlane.hpp 

# these will be installed
libmbcslam_la_include_HEADERS = 

libmbcslam_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = intx_on_sphere_test  intx_in_plane_test  spec_visu_test spherical_area_test \
         case1_test  intx_mpas
noinst_PROGRAMS =  cslam_par_test diffusion proj1 

if ENABLE_mbzoltan
noinst_PROGRAMS+=MPASpartitioner
MPASpartitioner_SOURCES = MPASpartitioner.cpp
endif 

check_PROGRAMS = $(TESTS) 
intx_on_sphere_test_SOURCES = intx_on_sphere_test.cpp
diffusion_SOURCES = diffusion.cpp
intx_in_plane_test_SOURCES = intx_in_plane_test.cpp
spec_visu_test_SOURCES = spec_visu_test.cpp
spherical_area_test_SOURCES = spherical_area_test.cpp
case1_test_SOURCES = case1_test.cpp
intx_mpas_SOURCES = intx_mpas.cpp
cslam_par_test_SOURCES = cslam_par_test.cpp
proj1_SOURCES = proj1.cpp

EXTRA_DIST  = lagrangeHomme.vtk  \
              eulerHomme.vtk \
              m1.vtk \
              m2.vtk 

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
           
