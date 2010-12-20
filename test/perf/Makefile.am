
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/src \
               -DIS_BUILDING_MB
            
LDADD = $(top_builddir)/src/libMOAB.la

check_PROGRAMS = perf seqperf adj_time perftool
noinst_PROGRAMS = point_in_elem 

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
adj_time_SOURCES = adj_time.cpp
point_in_elem_SOURCES = point_in_elem.cpp
perftool_SOURCES = perftool.cpp

if ENABLE_imesh
  LDADD += $(top_builddir)/itaps/imesh/libiMesh.la
  AM_CPPFLAGS += -I$(top_srcdir)/itaps/imesh -I$(top_srcdir)/itaps -I$(top_builddir)/itaps -I$(top_builddir)/itaps/imesh
  check_PROGRAMS += tstt_perf_binding
  tstt_perf_binding_DEPENDENCIES = $(top_builddir)/itaps/imesh/libiMesh.la
  tstt_perf_binding_SOURCES = tstt_perf_binding.cpp
endif


