AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src -I$(top_srcdir)/tools/mbcoupler -I$(top_builddir)/src/moab/point_locater/lotte

bin_PROGRAMS = mcnpmit
mcnpmit_SOURCES = mcnpmit.hpp mcnpmit.cpp main.cpp
mcnpmit_LDADD = $(top_builddir)/tools/mbcoupler/libmbcoupler.la $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la

