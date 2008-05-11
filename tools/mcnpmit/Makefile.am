AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = mcnpmit
mcnpmit_SOURCES = mcnpmit.hpp mcnpmit.cpp main.cpp
LDADD = $(top_builddir)/libMOAB.la
mcnpmit_DEPENDENCIES = $(top_builddir)/libMOAB.la
