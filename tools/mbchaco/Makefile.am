AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir)
bin_PROGRAMS = mbchaco
mbchaco_SOURCES = ComputePartition.hpp ComputePartition.cpp main.cpp
mbchaco_LDADD = $(LIBS)
mbchaco_DEPENDENCIES = $(top_builddir)/libMOAB.la
