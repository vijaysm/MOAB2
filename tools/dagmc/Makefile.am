AUTOMAKE_OPTIONS = foreign

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) 

lib_LTLIBRARIES = libdagmc.la
libdagmc_la_SOURCES = DagMC.cpp
libdagmc_la_includedir = $(includedir)
libdagmc_la_include_HEADERS = cgm2moab.hpp cubfile.h DagMC.hpp

if HAVE_CGM
INCLUDES += $(CGM_INCLUDES)
CPPFLAGS += $(CGM_CPPFLAGS) -DCGM $(MOAB_CGM_DEFINES)
CXXFLAGS += $(CGM_CXXFLAGS)

libdagmc_la_SOURCES += cgm2moab.cc

noinst_PROGRAMS = update_coords quads_to_tris
bin_PROGRAMS = cgm2moab
cgm2moab_SOURCES = main.cc cubfile.c
cgm2moab_LDFLAGS = $(CGM_LTFLAGS) $(CGM_LDLFAGS)
cgm2moab_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LIBS)
cgm2moab_DEPENDENCIES = libdagmc.la $(top_builddir)/libMOAB.la
dist_man1_MANS = cgm2moab.man
endif
                   
CGM_CPPFLAGS = @CGM_CPPFLAGS@
CGM_LDFLAGS = @CGM_LDFLAGS@
CGM_LTFLAGS = @CGM_LTFLAGS@
CGM_LIBS = @CGM_LIBS@

TESTS = test_geom 
check_PROGRAMS = $(TESTS) pt_vol_test ray_fire_test

quads_to_tris_SOURCES = quads_to_tris_driver.cpp quads_to_tris.cpp
quads_to_tris_LDADD = $(top_builddir)/libMOAB.la

test_geom_SOURCES = test_geom.cc
test_geom_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LTFLAGS)  $(CGM_LDFLAGS) $(CGM_LIBS)

update_coords_SOURCES = cub2h5m.cc quads_to_tris.cpp quads_to_tris.hpp
update_coords_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LTFLAGS) $(CGM_LDFLAGS) $(CGM_LIBS)

pt_vol_test_SOURCES = pt_vol_test.cc
pt_vol_test_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LTFLAGS)  $(CGM_LDFLAGS) $(CGM_LIBS)

ray_fire_test_SOURCES = ray_fire_test.cc
ray_fire_test_LDADD = libdagmc.la $(top_builddir)/libMOAB.la $(CGM_LTFLAGS)  $(CGM_LDFLAGS) $(CGM_LIBS)
