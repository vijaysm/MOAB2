MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src 

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/src/moab/point_locater/lotte 

lib_LTLIBRARIES = libmbcslam.la
libmbcslam_la_LIBADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/itaps/imesh/libiMesh.la

LDADD = libmbcslam.la 

libmbcslam_la_SOURCES = \
   Intx2Mesh.cpp Intx2Mesh.hpp Intx2MeshOnSphere.cpp Intx2MeshOnSphere.hpp \
     IntxUtils.cpp IntxUtils.hpp

# these will be installed
libmbcslam_la_include_HEADERS = 

libmbcslam_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = intx_on_sphere_test
check_PROGRAMS = $(TESTS) 
intx_on_sphere_test_SOURCES = intx_on_sphere_test.cpp

