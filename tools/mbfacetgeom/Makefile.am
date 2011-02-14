
MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB\
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src 

lib_LTLIBRARIES = libmbfacetgeom.la
libmbfacetgeom_la_LIBADD = $(top_builddir)/src/libMOAB.la 
libmbfacetgeom_la_includedir = $(includedir)
libmbfacetgeom_la_SOURCES = FBEngine.cpp SmoothCurve.cpp SmoothFace.cpp
libmbfacetgeom_la_include_HEADERS = FBEngine.hpp

TESTS = facet_test

check_PROGRAMS = $(TESTS)

LDADD = libmbfacetgeom.la

facet_test_SOURCES = mbfacet_test.cpp
facet_test_DEPENDENCIES = \
  $(LDADD)

