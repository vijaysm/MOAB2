AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) $(CGM_INCLUDES)
CPPFLAGS += $(TEMPLATE_DEFS_INCLUDED)

bin_PROGRAMS = cgm2moab
cgm2moab_SOURCES = main.cc cubfile.c

libcgm2moab_la_SOURCES = \
	cgm2moab.cc
                   
libMOAB_la_include_HEADERS = cgm2moab.hpp cubfile.h

lib_LTLIBRARIES = libcgm2moab.la

libMOAB_la_includedir = $(includedir)

LDADD = libcgm2moab.la $(top_builddir)/libMOAB.la $(CGM_LIBS_LINK)
cgm2moab_DEPENDENCIES = $(top_builddir)/libMOAB.la
dist_man1_MANS = cgm2moab.man
