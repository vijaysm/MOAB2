AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir)
bin_PROGRAMS = mbconvert
mbconvert_SOURCES = convert.cpp
LDADD = $(top_builddir)/libMOAB.la
mbconvert_DEPENDENCIES = $(top_builddir)/libMOAB.la
