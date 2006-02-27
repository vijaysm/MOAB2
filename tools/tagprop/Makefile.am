AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir)
bin_PROGRAMS = mbtagprop
mbtagprop_SOURCES = parse.cpp parse.hpp propagate_tags.cpp
LDADD = $(top_builddir)/libMOAB.la
mbtagprop_DEPENDENCIES = $(top_builddir)/libMOAB.la
