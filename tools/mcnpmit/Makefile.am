AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)

LDADD = $(top_builddir)/libMOAB.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la

INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/tools/mbcoupler
bin_PROGRAMS = mcnpmit
mcnpmit_SOURCES = mcnpmit.hpp mcnpmit.cpp main.cpp
mcnpmit_DEPENDENCIES = $(top_builddir)/libMOAB.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la
