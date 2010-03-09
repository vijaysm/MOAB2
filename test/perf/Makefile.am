
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/moab \
               -I$(top_builddir)/src/moab
            
LDADD = $(top_builddir)/src/libMOAB.la

check_PROGRAMS = perf seqperf adj_time
noinst_PROGRAMS = point_in_elem

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
adj_time_SOURCES = adj_time.cpp
point_in_elem_SOURCES = point_in_elem.cpp


