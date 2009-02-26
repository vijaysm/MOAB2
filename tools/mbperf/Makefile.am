AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB

if USE_BABEL
  babel_includes = -I$(BABEL_DIR)/include
  babel_progs = mbperf_SIDL
else
  babel_includes = 
  babel_progs = 
endif

if ENABLE_imesh
  imesh_progs = mbperf_iMesh
else
  imesh_progs = 
endif 


INCLUDES += -I$(top_srcdir) \
            -I$(top_builddir) \
            -I$(top_srcdir)/tools/iMesh \
            -I$(top_srcdir)/tools/iMesh/SIDL/mserver \
            -I../iMesh/SIDL/mserver \
            -I../iMesh \
            $(babel_includes)

bin_PROGRAMS = mbperf $(babel_progs) $(imesh_progs)

mbperf_SOURCES = mbperf.cpp
mbperf_SIDL_SOURCES = mbperf_SIDL.cpp
mbperf_iMesh_SOURCES = mbperf_iMesh.cpp

LDADD = $(top_builddir)/libMOAB.la
mbperf_iMesh_LDADD = $(LDADD) ../iMesh/libiMesh.la 
mbperf_SIDL_LDADD = $(mbperf_iMesh_LDADD) ../iMesh/SIDL/mserver/libiMeshserver.la


