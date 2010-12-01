AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_builddir)/src 

lib_LTLIBRARIES = libdagmc.la
libdagmc_la_SOURCES = DagMC.cpp
libdagmc_la_includedir = $(includedir)
libdagmc_la_include_HEADERS = cubfile.h DagMC.hpp

if HAVE_CGM
  AM_CPPFLAGS += $(CGM_INCLUDES) $(CGM_CPPFLAGS) -DCGM $(MOAB_CGM_DEFINES)
  AM_CXXFLAGS += $(CGM_CXXFLAGS)

if NETCDF_FILE
  netcdf_progs = update_coords 
else
  netcdf_progs = 
endif

  cgm2moabbin = cgm2moab
  cgm2moabman = cgm2moab.man
else
  cgm2moabbin=
  cgm2moabman=
endif

cgm2moab_SOURCES = main.cc cubfile.c
cgm2moab_LDFLAGS = $(CGM_LTFLAGS) $(CGM_LDLFAGS)
cgm2moab_LDADD = $(LDADD) 
cgm2moab_DEPENDENCIES = libdagmc.la $(top_builddir)/src/libMOAB.la

bin_PROGRAMS = $(cgm2moabbin)
noinst_PROGRAMS = quads_to_tris $(netcdf_progs) dagmc_preproc 
man1_MANS = $(cgm2moabman)
EXTRA_DIST = cgm2moab.man

CGM_CPPFLAGS = @CGM_CPPFLAGS@
CGM_LDFLAGS = @CGM_LDFLAGS@
CGM_LTFLAGS = @CGM_LTFLAGS@
CGM_LIBS = @CGM_LIBS@

LDADD = libdagmc.la $(top_builddir)/src/libMOAB.la  $(CGM_LDFLAGS) $(CGM_LIBS)

TESTS = test_geom 
check_PROGRAMS = $(TESTS) pt_vol_test ray_fire_test

quads_to_tris_SOURCES = quads_to_tris_driver.cpp quads_to_tris.cpp
quads_to_tris_LDADD = $(top_builddir)/src/libMOAB.la

test_geom_SOURCES = test_geom.cc

update_coords_SOURCES = cub2h5m.cc quads_to_tris.cpp quads_to_tris.hpp

pt_vol_test_SOURCES = pt_vol_test.cc

ray_fire_test_SOURCES = ray_fire_test.cc

dagmc_preproc_SOURCES = dagmc_preproc.cpp ../program_opt.cpp
