LDADD = $(top_builddir)/src/libMOAB.la

if ENABLE_imesh
  imesh_progs = mbperf_iMesh
else
  imesh_progs = 
endif 


AM_CPPFLAGS += -DIS_BUILDING_MB \
            -I$(top_srcdir)/src \
            -I$(top_builddir)/src \
            -I$(top_srcdir)/itaps \
            -I$(top_srcdir)/itaps/imesh \
            -I$(top_builddir)/itaps \
            -I$(top_builddir)/itaps/imesh \
            -I$(top_srcdir)/src/moab \
            -I$(top_builddir)/src/moab

bin_PROGRAMS = mbperf $(imesh_progs)

mbperf_SOURCES = mbperf.cpp
mbperf_iMesh_SOURCES = mbperf_iMesh.cpp

mbperf_iMesh_LDADD = $(LDADD) $(top_builddir)/itaps/imesh/libiMesh.la 


