AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
bin_PROGRAMS = mbconvert
mbconvert_SOURCES = convert.cpp
LDADD = $(top_builddir)/libMOAB.la
mbconvert_DEPENDENCIES = $(top_builddir)/libMOAB.la
dist_man1_MANS = mbconvert.man
