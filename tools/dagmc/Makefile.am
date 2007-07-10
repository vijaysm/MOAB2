AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) 

lib_LTLIBRARIES = libdagmc.la
libdagmc_la_SOURCES = DagMC.cpp
libdagmc_la_includedir = $(includedir)
libdagmc_la_include_HEADERS = cgm2moab.hpp cubfile.h DagMC.hpp

if HAVE_CGM
INCLUDES += $(CGM_INCLUDES)
CPPFLAGS += $(CGM_CPPFLAGS) -DCGM
CXXFLAGS += $(CGM_CXXFLAGS)

libdagmc_la_SOURCES += cgm2moab.cc

bin_PROGRAMS = cgm2moab
cgm2moab_SOURCES = main.cc cubfile.c
cgm2moab_LDFLAGS = $(CGM_LIBS_LTFLAGS)
cgm2moab_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LIBS_LINK)
cgm2moab_DEPENDENCIES = libdagmc.la $(top_builddir)/libMOAB.la
dist_man1_MANS = cgm2moab.man
endif
                   

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = dagmc.make

TESTS = test_geom
check_PROGRAMS = $(TESTS)
test_geom_SOURCES = test_geom.cc
test_geom_LDADD = libdagmc.la $(top_builddir)/libMOAB.la


