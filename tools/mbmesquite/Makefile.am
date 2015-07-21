MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src 

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/tools \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/test

lib_LTLIBRARIES = libmbmesquite.la
LDADD = libmbmesquite.la  $(top_builddir)/src/libMOAB.la $(MESQUITE_LIBS)
AM_CPPFLAGS += $(MESQUITE_INC_FLAGS)
AM_LDFLAGS += $(MESQUITE_LIB_FLAGS)
TESTS = imesh_mesquite igeom_mesquite


if ENABLE_imesh
 LDADD+= $(top_builddir)/itaps/imesh/libiMesh.la
endif

libmbmesquite_la_SOURCES = \
   MsqIMesh.cpp MsqIMesh.hpp MsqIBase.cpp MsqIBase.hpp

if HAVE_CGM
  AM_CPPFLAGS += $(CGM_CPPFLAGS)
  LDADD+=  ${CGM_DIR}/lib/libiGeom.la ${CGM_DIR}/lib/libcgm.la
  libmbmesquite_la_SOURCES += MsqIGeom.cpp MsqIGeom.hpp
endif

if PARALLEL
    libmbmesquite_la_SOURCES += MsqIMeshP.cpp MsqIMeshP.hpp
endif


# these will be installed
libmbmesquite_la_include_HEADERS =

libmbmesquite_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

check_PROGRAMS = $(TESTS)

imesh_mesquite_SOURCES = imesh_mesquite.cpp
igeom_mesquite_SOURCES = igeom_mesquite.cpp

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = smart-laplacian-result.vtk \
                   feasible-newton-result.vtk \
                   original.vtk \
                   test-suite.log
