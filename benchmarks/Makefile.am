AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB

LDADD = $(top_builddir)/libMOAB.la

INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/parallel \
            -I$(top_builddir)
bin_PROGRAMS = mbparead
mbparead_SOURCES = mbparead.cpp 
mbparead_DEPENDENCIES = $(top_builddir)/libMOAB.la


