AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir)
bin_PROGRAMS = mbsize
mbsize_SOURCES = measure.cpp measure.hpp size.cpp
LDADD = $(top_builddir)/libMOAB.la
mbsize_DEPENDENCIES = $(top_builddir)/libMOAB.la
