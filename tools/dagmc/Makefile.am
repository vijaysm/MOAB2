AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) $(CGM_INCLUDES)
CPPFLAGS += $(CGM_CPPFLAGS)
CXXFLAGS += $(CGM_CXXFLAGS)

bin_PROGRAMS = cgm2moab
cgm2moab_SOURCES = main.cc cubfile.c
cgm2moab_LDFLAGS = $(CGM_LIBS_LTFLAGS)

libdagmc_la_SOURCES = \
	cgm2moab.cc DagMC.cpp
                   
libMOAB_la_include_HEADERS = cgm2moab.hpp cubfile.h DagMC.hpp

lib_LTLIBRARIES = libdagmc.la

libMOAB_la_includedir = $(includedir)

LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LIBS_LINK)
cgm2moab_DEPENDENCIES = libdagmc.la $(top_builddir)/libMOAB.la
dist_man1_MANS = cgm2moab.man

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = dagmc.make

