AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)

@CGM_CONFIG_OPTIONS@
LDADD = $(top_builddir)/libMOAB.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la
if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif

INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/tools/mbcoupler
bin_PROGRAMS = mcnpmit
mcnpmit_SOURCES = mcnpmit.hpp mcnpmit.cpp main.cpp
mcnpmit_DEPENDENCIES = $(top_builddir)/libMOAB.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la
