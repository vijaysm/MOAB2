AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_builddir)/src $(CHACO_INCLUDES) \
               -I$(top_srcdir)/src/moab -I$(top_builddir)/src/moab
bin_PROGRAMS = mbchaco
mbchaco_SOURCES = ComputePartition.hpp ComputePartition.cpp main.cpp
LDADD = $(top_builddir)/src/libMOAB.la $(CHACO_LIBS)

mbchaco_DEPENDENCIES = $(top_builddir)/src/libMOAB.la
