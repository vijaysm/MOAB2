AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir) -I$(top_builddir)
bin_PROGRAMS = mbsurfplot
mbsurfplot_SOURCES = surfplot.cpp
LDADD = $(top_builddir)/libMOAB.la
mbsurfplot_DEPENDENCIES = $(top_builddir)/libMOAB.la
