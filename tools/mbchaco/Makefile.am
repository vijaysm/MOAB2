AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) $(CHACO_INCLUDES)
bin_PROGRAMS = mbchaco
mbchaco_SOURCES = ComputePartition.hpp ComputePartition.cpp main.cpp
LDADD = $(top_builddir)/libMOAB.la $(CHACO_LIBS)

if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif

mbchaco_DEPENDENCIES = $(top_builddir)/libMOAB.la
