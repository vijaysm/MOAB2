AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = mbgsets
mbgsets_SOURCES = gsets.cc
LDADD = $(top_builddir)/libMOAB.la
mbgsets_DEPENDENCIES = $(top_builddir)/libMOAB.la
