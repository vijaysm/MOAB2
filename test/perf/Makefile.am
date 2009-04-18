
@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/mhdf/include \
            -I$(top_builddir) \
            -I$(top_srcdir)/tools/iMesh \
            -I$(top_builddir)/tools/iMesh
            
LDADD = $(top_builddir)/libMOAB.la

if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif

check_PROGRAMS = perf seqperf adj_time
noinst_PROGRAMS = point_in_elem

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
adj_time_SOURCES = adj_time.cpp
point_in_elem_SOURCES = point_in_elem.cpp


