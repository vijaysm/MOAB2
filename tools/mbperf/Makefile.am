AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB

INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = mbperf
mbperf_SOURCES = mbperf.cpp
LDADD = $(top_builddir)/libMOAB.la
mbperf_DEPENDENCIES = $(top_builddir)/libMOAB.la

if USE_BABEL
  INCLUDES += -I$(top_srcdir)/tools/iMesh/SIDL/mserver \
              -I../iMesh/SIDL/mserver \
              -I$(BABEL_DIR)/include
  bin_PROGRAMS += mbperf_SIDL
  mbperf_SIDL_SOURCES = mbperf_SIDL.cpp
  LDADD += $(top_builddir)/tools/iMesh/SIDL/mserver/libiMeshserver.la
endif

if ENABLE_imesh
  INCLUDES += -I$(top_srcdir)/tools/iMesh/
  bin_PROGRAMS += mbperf_iMesh
  LDADD += $(top_builddir)/tools/iMesh/libiMesh.la 
  mbperf_iMesh_SOURCES = mbperf_iMesh.cpp
endif 

