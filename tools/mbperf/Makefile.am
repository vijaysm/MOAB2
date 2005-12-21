AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB

INCLUDES += -I$(top_builddir)
bin_PROGRAMS = mbperf
mbperf_SOURCES = mbperf.cpp
LDADD = $(top_builddir)/libMOAB.la
mbperf_DEPENDENCIES = $(top_builddir)/libMOAB.la
