AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) $(CGM_INCLUDES)
bin_PROGRAMS = cgm2moab
cgm2moab_SOURCES = cgm2moab.cc \
                   cubfile.c \
                   cubfile.h \
                   geom_file_type.c \
                   geom_file_type.h
                   
LDADD = $(top_builddir)/libMOAB.la $(CGM_LIBS_LINK)
cgm2moab_DEPENDENCIES = $(top_builddir)/libMOAB.la
dist_man1_MANS = cgm2moab.man
