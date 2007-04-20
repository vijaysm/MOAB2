AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = hexmodops
hexmodops_SOURCES = makeops.cpp
LDADD = $(top_builddir)/libMOAB.la
hexmodops_DEPENDENCIES = $(top_builddir)/libMOAB.la
