AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = spheredecomp
spheredecomp_SOURCES = SphereDecomp.hpp SphereDecomp.cpp main.cpp
LDADD = $(top_builddir)/libMOAB.la
spheredecomp_DEPENDENCIES = $(top_builddir)/libMOAB.la
