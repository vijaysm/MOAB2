DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/mhdf/include \
            -I$(top_builddir) \
            -I$(top_srcdir)/tools/iMesh \
            -I$(top_builddir)/tools/iMesh
            
LDADD = $(top_builddir)/libMOAB.la

check_PROGRAMS = perf seqperf adj_time

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
adj_time_SOURCES = adj_time.cpp

