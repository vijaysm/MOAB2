AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir)
bin_PROGRAMS = hexmodops
hexmodops_SOURCES = makeops.cpp
LDADD = $(top_builddir)/libMOAB.la
hexmodops_DEPENDENCIES = $(top_builddir)/libMOAB.la
