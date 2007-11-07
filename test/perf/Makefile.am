DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/mhdf/include \
            -I$(top_builddir) \
            -I$(top_srcdir)/tools/iMesh \
            -I$(top_builddir)/tools/iMesh
            
LDADD = $(top_builddir)/libMOAB.la

check_PROGRAMS = perf seqperf #cubit_perf tstt_perf tstt_perf_binding 

perf_SOURCES = perf.cpp
seqperf_SOURCES = seqperf.cpp
#cubit_perf_SOURCES = cubit_perf.cpp
#tstt_perf_SOURCES = tstt_perf.cpp
#tstt_perf_binding_SOURCES = tstt_perf_binding.cpp

