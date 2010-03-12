LDADD = $(top_builddir)/src/libMOAB.la \
        $(top_builddir)/tools/mbcoupler/libmbcoupler.la

AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_builddir)/src -I$(top_srcdir)/tools/mbcoupler

bin_PROGRAMS = mcnpmit
mcnpmit_SOURCES = mcnpmit.hpp mcnpmit.cpp main.cpp
mcnpmit_DEPENDENCIES = ../../src/libMOAB.la ../mbcoupler/libmbcoupler.la
