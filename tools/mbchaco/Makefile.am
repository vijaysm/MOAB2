AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir)
bin_PROGRAMS = mbchaco
mbchaco_SOURCES = ComputePartition.hpp ComputePartition.cpp main.cpp
LDADD = $(top_builddir)/../../lib/libMOAB.la
mbchaco_DEPENDENCIES = $(top_builddir)/../../lib/libMOAB.la
